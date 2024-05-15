import streamlit as st
import pandas as pd
import os

def runUI():
    if not st.session_state["queue"]:
        st.session_state["queue"] = True

    def get_job_example():
        st.session_state["job_input"] = "aGtWEUgLgMEvctFh"

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

        predictions = os.path.join(st.session_state["job_path"], "model_predictions.csv")
        df_results = pd.read_csv(predictions)

        df_results["Probability"] = df_results.apply(lambda x: max(x[["Cis-reg", "coding", "rRNA", "sRNA", "tRNA", "unknown"]])*100, axis=1)

        df_results = df_results[["nameseq", "prediction", "Probability"]]
        
        df_results.columns = ["Name", "Prediction", "Probability"]

        st.dataframe(df_results,
                    column_config = {"Probability": st.column_config.ProgressColumn(
                        help="Prediction probability",
                        format="%.2f%%",
                        min_value=0,
                        max_value=100
                    )},
                    hide_index=True, use_container_width=True)