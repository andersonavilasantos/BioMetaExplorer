import streamlit as st
import polars as pl
from math import ceil
import plotly.express as px

@st.cache_data(show_spinner=False)
def fetch_data():
    df_pred = pl.read_csv("../Classification/data/predict/predictions_complete.tsv", separator="\t")

    return df_pred.to_pandas()

@st.cache_data(show_spinner=False)
def convert_df(df):
    return df.to_csv(sep="\t").encode('utf-8')

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

            page["Probability"] = page.apply(lambda x: max(x[["Cis-reg", "coding", "rRNA", "sRNA", "tRNA", "unknown"]]), axis=1)

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
                                        help="Prediction probability",
                                        format="%.2f",
                                        min_value=0,
                                        max_value=1
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

        # classes = ["Cis-reg", "coding", "rRNA", "sRNA", "tRNA", "unknown"]

        # fig = px.pie(values=selected_rows[classes], names=classes)

        # st.plotly_chart(fig, use_container_width=True, help="Taxonomy")