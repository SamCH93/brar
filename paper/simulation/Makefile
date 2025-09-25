all: simulation-dashboard.html

DATA := $(wildcard  results/*.rds)

sim-summaries.rds: simulation-summaries.R $(DATA)
	Rscript simulation-summaries.R

simulation-dashboard.html: simulation-dashboard.Qmd sim-summaries.rds
	quarto render simulation-dashboard.Qmd --output-dir dashboard
