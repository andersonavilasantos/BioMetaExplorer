import streamlit as st
import polars as pl
import pandas as pd
from math import ceil
import numpy as np
from plotly.subplots import make_subplots
import plotly.express as px
import plotly.graph_objects as go
import matplotlib.pyplot as plt
from pyfamsa import Aligner, Sequence
import utils
import RNA

@st.cache_data(show_spinner=False)
def fetch_data():
    df_pred = pl.read_csv("../Classification/data/predict/predictions_complete.tsv", separator="\t")

    return df_pred.to_pandas()

@st.cache_data(show_spinner=False)
def convert_df(df):
    return df.to_csv(sep="\t").encode('utf-8')

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
    # df_sequences = pd.DataFrame()

    # for seq_class in seqs:
    #     df_seq = seqs[seq_class].df.copy()
    #     df_seq['list_name'] = seq_class + ' - ' + df_seq['name']
    #     df_sequences = pd.concat([df_sequences, df_seq]).reset_index(drop=True)

    # seq_select = st.multiselect("Select sequences to view alignment:", df_sequences["nameseq"])

    chars = ['A', 'C', 'G', 'T']

    # if seq_select:
    with st.spinner('Loading...'):
        # nseqs_selected = len(seq_select)

        # df_select = df_sequences[df_sequences["nameseq"].isin(seq_select)]

        sequences = [Sequence(name.encode(), 
                    df_sequences.loc[df_sequences["nameseq"] == name]["Sequence"].item().encode()) 
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
        seq_select = st.selectbox("Select sequence to predict structure:", 
                                options=df_sequences["nameseq"], index=None,
                                placeholder="Choose a sequence")

    if seq_select:
        seq = df_sequences[df_sequences["nameseq"] == seq_select].reset_index(drop=True)["Sequence"][0]

        ss, _ = RNA.fold(seq)

        # home_dir = os.path.expanduser('~')
        # dir_path = os.path.join(home_dir, '.biotukey')
        
        RNA.svg_rna_plot(seq, ss, "rna_plot.svg")

        with col2:
            st.markdown("**Dot-bracket notation:**")
            st.markdown(ss)
            st.markdown("**Secondary structure:**")
            st.image("rna_plot.svg", use_column_width = 'always')

def runUI():
    data = fetch_data()

    data_size = len(data)

    # download = st.columns([4, 1])

    # with download[1]:
    #     st.download_button(
    #         label="Download complete data as TSV",
    #         data=convert_df(data),
    #         file_name="data.tsv",
    #         mime='text/tsv',
    #         use_container_width=True
    #     )

    table = st.container()

    page_size = 100

    menu = st.columns([6, 1])

    total_pages = ceil(data_size/page_size)

    with menu[1]:

        page_number = st.number_input(
            label="Page",
            label_visibility="collapsed",
            min_value=1,
            max_value=total_pages,
            step=1,
        )

    with menu[0]:
        st.markdown(f"Page **{page_number}** of **{total_pages}** ")

    current_start = (page_number - 1) * page_size

    with table:
        with st.spinner("Loading..."):
            page = data.iloc[current_start:current_start + page_size, :].copy()

            page["Probability"] = page.apply(lambda x: max(x[["Cis-reg", "coding", "rRNA", "sRNA", "tRNA", "unknown"]])*100, axis=1)

            show_columns = ["mag", "GTDB-tk_domain", "GTDB-tk_phylum", "prediction"]

            col1, col2 = st.columns([2, 1])

            with col1:
                page.insert(0, "View", False)

                edited_df = st.data_editor(
                    page,
                    hide_index=True,
                    height=500,
                    column_config = {"View": st.column_config.CheckboxColumn(required=True),
                                    "Probability": st.column_config.ProgressColumn(
                                        help="Probability of predicted class",
                                        format="%.2f%%",
                                        min_value=0,
                                        max_value=100
                                    ),},
                    column_order=["View"] + show_columns + ["Probability"],
                    disabled=show_columns,
                    use_container_width=True
                )

                selected_rows = edited_df[edited_df["View"]].drop("View", axis=1).reset_index(drop=True)
                
            with col2:
                fig = px.sunburst(
                    page,
                    path=["GTDB-tk_domain", "GTDB-tk_phylum", "prediction"]
                    # parents="GTDB-tk_domain",
                    # names="GTDB-tk_phylum",
                    
                    # values='nameseq',
                )
                
                st.plotly_chart(fig, use_container_width=True, help="Taxonomy")
    
    if not selected_rows.empty:
        st.divider()

        st.dataframe(selected_rows, hide_index=True)

        tab1, tab2 = st.tabs(["Sequence Alignment", "Secondary Structure"])

        with tab1:
            seq_alignment(selected_rows)

        with tab2:
            structure_visualization(selected_rows)

        # classes = ["Cis-reg", "coding", "rRNA", "sRNA", "tRNA", "unknown"]

        # fig = px.pie(values=selected_rows[classes], names=classes)

        # st.plotly_chart(fig, use_container_width=True, help="Taxonomy")