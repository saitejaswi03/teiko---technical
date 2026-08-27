.PHONY: setup pipeline dashboard clean

# Installs all dependencies needed to run the pipeline and dashboard.
setup:
	pip install --upgrade pip
	pip install -r requirements.txt

# Runs the full pipeline end to end with no manual intervention:
# initializes + loads the database (Part 1), then generates all Part 2-4
# output tables and figures into output/.
pipeline:
	python pipeline.py

# Starts the local server for the interactive Streamlit dashboard.
dashboard:
	streamlit run dashboard/app.py

# Convenience target (not required by the assignment): removes generated
# artifacts so `make pipeline` can be re-run from a clean state.
clean:
	rm -f cell_counts.db
	rm -rf output/*.csv output/*.png output/*.txt
