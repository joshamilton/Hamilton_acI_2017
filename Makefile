mansucript: manuscript/figureCaptions.md manuscript/manuscript.md manuscript/som.md
		pandoc -s -S manuscript/manuscript.md -o manuscript/Hamilton_acI_2016_MS.docx --reference-docx=manuscript/template.docx --filter pandoc-citeproc --bibliography manuscript/reverseEcology.bib --csl=manuscript/the-isme-journal.csl
		pandoc -s -S manuscript/figureCaptions.md -o manuscript/Hamilton_acI_2016_Captions.docx --reference-docx=manuscript/template.docx --filter pandoc-citeproc --bibliography manuscript/reverseEcology.bib --csl=manuscript/the-isme-journal.csl
		pandoc -s -S manuscript/som.md -o manuscript/Hamilton_acI_2016_SOM.docx --reference-docx=manuscript/template.docx --filter pandoc-citeproc --bibliography manuscript/reverseEcology.bib --csl=manuscript/the-isme-journal.csl
