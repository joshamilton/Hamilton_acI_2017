mansucript: manuscript/figureCaptions.md manuscript/manuscript.md
		pandoc -s -S manuscript/figureCaptions.md -o manuscript/figureCaptions.docx --reference-docx=manuscript/template.docx --filter pandoc-citeproc --bibliography manuscript/reverseEcology.bib -csl=manuscript/the-isme-journal.csl
		pandoc -s -S manuscript/manuscript.md -o manuscript/Hamilton_ISMEJ_2016.docx --reference-docx=manuscript/template.docx --filter pandoc-citeproc --bibliography manuscript/reverseEcology.bib -csl=manuscript/the-isme-journal.csl
