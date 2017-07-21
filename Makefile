mansucript: manuscript/manuscript.md manuscript/som.md
		pandoc -s -S manuscript/manuscript.md -o manuscript/Hamilton_acI_2017_MS.docx --reference-docx=manuscript/template.docx --filter pandoc-citeproc --bibliography manuscript/reverseEcology.bib --csl=manuscript/american-society-for-microbiology.csl
		pandoc -s -S manuscript/som.md -o manuscript/Hamilton_acI_2017_SOM.docx --reference-docx=manuscript/template.docx --filter pandoc-citeproc --bibliography manuscript/reverseEcology.bib --csl=manuscript/american-society-for-microbiology.csl
