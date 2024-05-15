import plotly.express as px
import plotly.graph_objects as go
import streamlit as st
from pyfamsa import Aligner, Sequence
from plotly.subplots import make_subplots
import numpy as np
import RNA
import pandas as pd
import os
import time
import utils
import tempfile

def information_content(col, align_length):
    rel_freqs = [b/sum(col[1:]) for b in col[1:]]
    s = len(rel_freqs)

    entropy = -sum([freq * np.log2(freq) if freq != 0 else 0 for freq in rel_freqs])

    sample_correction = (1 / np.log(2)) * ((s - 1)/ (2 * align_length))
    R_i = np.log2(s) - (entropy + sample_correction)

    return R_i

def relative_information(col, freq, ic):
    rel_freq = freq/sum(col[1:])

    height = rel_freq * ic

    return height

def seq_alignment(df_sequences):
    chars = ['A', 'C', 'G', 'T']

    # if seq_select:
    with st.spinner('Loading...'):
        # nseqs_selected = len(seq_select)

        # df_select = df_sequences[df_sequences["nameseq"].isin(seq_select)]

        sequences = [Sequence(name.encode(), 
                    df_sequences.loc[df_sequences["nameseq"] == name]["sequence"].item().encode()) 
                    for name in df_sequences["nameseq"]]

        aligner = Aligner(guide_tree="upgma")
        msa = aligner.align(sequences)
        
        msa_dict = {sequence.id.decode():sequence.sequence.decode() for sequence in msa}
        
        char_colors = {char: i / len(chars) for i, char in enumerate(chars)}
        char_colors['-'] = 1

        msa_seqs = [[*msa_dict[id]] for id in msa_dict]

        seqs_num = [[char_colors[N] for N in msa_seq] for msa_seq in msa_seqs]

        colorscale = [[char_colors[ch_color], utils.get_colors(len(char_colors))[i]] for i, ch_color in enumerate(char_colors)]
        colorscale[-1] = [1, 'white']
        
        fig = make_subplots(rows=2, shared_xaxes=True, vertical_spacing=0.065, row_heights=[0.2, 0.8])

        np_seqs = np.array(msa_seqs)

        len_seqs = len(np_seqs[0])

        df_dict = {'Position': list(range(len_seqs))}

        for char in chars:
            df_dict[char] = [list(np_seqs[:,i]).count(char) for i in range(len_seqs)]

        df = pd.DataFrame(df_dict)

        colors = {char: utils.get_colors(len(chars))[i] for i, char in enumerate(chars)}

        data = []
        
        for i in range(df.shape[0]):
            ordered_columns = df.columns[1:][np.argsort(df.iloc[i, 1:].values)]

            info_col = information_content(df.iloc[i], df.shape[0])
            for _, column in enumerate(ordered_columns):
                data.append(go.Bar(x=[df['Position'][i]],
                                y=[relative_information(df.iloc[i], df[column][i], info_col)],
                                textfont=dict(
                                    family="Trebuchet MS",
                                    size = 10,
                                ),
                                cliponaxis = False,
                                marker=dict(color=colors[column]),
                                name = column,
                                textposition="outside",
                                legendgroup=column,
                                showlegend=i == 0,
                                hovertemplate = "<b>Nucleotide: </b>" + column + " "
                                        "<br><b>Information content:</b> %{y}" +
                                        "<br><b>Position:</b> %{x} <extra></extra>")) 

        fig_hm = px.imshow(seqs_num, aspect = 'auto')
        fig_hm.update_traces(text=msa_seqs,
                        hovertemplate = "<b>Sequence name:</b> %{y}" +
                                        "<br> <b>Nucleotide:</b> %{text}" +
                                        "<br> <b>Position:</b> %{x} <extra></extra>",
                        
                        textfont=dict(
                                family="Trebuchet MS",
                                size = 6,
                            ), texttemplate="<b>%{text}</b>")
        
        for bar in data:
            fig.add_trace(bar, row=1, col=1)

        fig.add_trace(fig_hm['data'][0], row = 2, col = 1)

        fig.update_coloraxes(colorscale = colorscale, cmin = 0, cmax = 1, showscale=False)
        fig.update_layout(height = 800, barmode = 'stack', dragmode='pan',
                        legend_title_text='Nucleotides',
                        yaxis2 = dict(tickmode = 'array',
                            tickvals = list(range(len(df_sequences))),
                            ticktext = list(msa_dict.keys())),
                        xaxis2=dict(
                        rangeslider=dict(visible=True), range = (-0.5, 50), type="linear",))

        st.plotly_chart(fig, use_container_width=True)

def structure_visualization(df_sequences):
    col1, col2 = st.columns(2)

    with col1:
        with st.form("structure_submit"):
            seq_select = st.selectbox("Select sequence to predict structure:", 
                        options=df_sequences["nameseq"], index=None,
                        placeholder="Choose a sequence")

            submitted = st.form_submit_button("Submit")

    if submitted:
        seq = df_sequences[df_sequences["nameseq"] == seq_select].reset_index(drop=True)["sequence"][0]

        with tempfile.TemporaryDirectory() as temp_dir:

            ss, _ = RNA.fold(seq)
            tmp_plot = os.path.join(temp_dir, f"{int(time.time())}.svg")
            RNA.svg_rna_plot(seq, ss, tmp_plot)

            with col2:
                st.markdown("**Dot-bracket notation:**")
                st.markdown(ss)
                st.markdown("**Structure visualization:**")
                st.image(tmp_plot, use_column_width = "always")

def runUI():
    if not st.session_state["queue"]:
        st.session_state["queue"] = True

    def get_job_example():
        st.session_state["job_input"] = "aqx2sFL35TmQdIcD"

    with st.container(border=True):
        col1, col2 = st.columns([9, 1])

        with col2:
            example_job = st.button("Example", use_container_width=True, on_click=get_job_example)

        with st.form("jobs_submit", border=False):
            job_id = st.text_input("Enter Job ID", key="job_input")

            submitted = st.form_submit_button("Submit", use_container_width=True,  type="primary")

    predict_path = "jobs"

    if submitted:
        if job_id:
            job_path = os.path.join(predict_path, job_id)
            if os.path.exists(job_path):
                predictions = os.path.join(job_path, "model_predictions.csv")
                if os.path.exists(predictions):
                    st.session_state["job_path"] = job_path
                else:
                    if "job_path" in st.session_state:
                        del st.session_state["job_path"]
                    st.info("Job is still in progress. Come back later.")
            else:
                if "job_path" in st.session_state:
                    del st.session_state["job_path"]
                st.error("Job does not exist!")

    if "job_path" in st.session_state:
        st.success("Job was completed with the following results.")

        with st.spinner("Loading..."):
            stats = os.path.join(st.session_state["job_path"], "stats.csv")
            df_stats = pd.read_csv(stats)

            st.dataframe(df_stats, hide_index=True, use_container_width=True)

            predictions = os.path.join(st.session_state["job_path"], "model_predictions.csv")
            df_results = pd.read_csv(predictions)

            df_results["probability"] = df_results.apply(lambda x: max(x[["Cis-reg", "coding", "rRNA", "sRNA", "tRNA", "unknown"]])*100, axis=1)

            show_columns = ["nameseq", "prediction"]

            col1, col2 = st.columns([3, 2])

            with col1:
                df_results.insert(0, "View", False)

                edited_df = st.data_editor(
                    df_results,
                    hide_index=True,
                    height=500,
                    column_config = {"View": st.column_config.CheckboxColumn(required=True),
                                    "probability": st.column_config.ProgressColumn(
                                        help="Prediction probability",
                                        format="%.2f%%",
                                        min_value=0,
                                        max_value=100
                                    ),},
                    column_order=["View"] + show_columns + ["probability"],
                    disabled=show_columns,
                    use_container_width=True
                )

                selected_rows = edited_df[edited_df["View"]].drop("View", axis=1).reset_index(drop=True)

            with col2:
                page = df_results.groupby(["prediction"])["prediction"].count().reset_index(name="count")

                fig =  px.pie(page, values="count", names="prediction")

                st.plotly_chart(fig, use_container_width=True, help="Taxonomy")

        if not selected_rows.empty:
            st.divider()

            st.dataframe(selected_rows.drop(columns=["probability"]), hide_index=True, use_container_width=True)

            tab1, tab2 = st.tabs(["Sequence Alignment", "Secondary Structure"])

            with tab1:
                seq_alignment(selected_rows)

            with tab2:
                structure_visualization(selected_rows)