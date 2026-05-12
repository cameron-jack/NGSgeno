#!/bin/bash

# you must have google-chrome and python 3.12+ installed
# We need a script to convert web links to #targets so that you can click through the PDF

python -m http.server 8000 &
for f in complete_NGSG_doc.html
do
    fx=${f/.html/.pdf}
    google-chrome --headless --print-to-pdf="${fx}" --no-pdf-header-footer http://localhost:8000/"${f}"
done

# after this you can combine the individual pdfs online with Adobe online for free

