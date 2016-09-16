figureCaptions.docx: manuscript/figureCaptions.md
		pandoc -s -S manuscript/figureCaptions.md -o manuscript/figureCaptions.docx --reference-docx=manuscript/template.docx --filter pandoc-citeproc --bibliography manuscript/reverseEcology.bib -csl=manuscript/the-isme-journal.csl
