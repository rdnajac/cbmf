# Makefile for building a LaTeX document and cleaning up auxiliary files

# Define the name of your LaTeX document (without the .tex extension)
DOC = your_document_name

# Define the LaTeX compiler
LATEXMK = latexmk

# Define the flags for latexmk
LATEXMKFLAGS = -pdf -silent

# Define the clean-up command
CLEANUP = rm -f *.aux *.log *.fls *.fdb_latexmk *.synctex.gz *.out

# Default target to build the PDF
all: $(DOC).pdf

# Rule to build the PDF using latexmk
$(DOC).pdf: $(DOC).tex
	$(LATEXMK) $(LATEXMKFLAGS) $(DOC).tex

# Clean up auxiliary files
clean:
	$(CLEANUP)

