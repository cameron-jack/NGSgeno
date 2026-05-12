#!/usr/bin/env python3
"""
This scripts reads the html2doc.cfg file to determine which HTML files to combine into a single html file, and then removes the headers and footers from the middle files.
generate_pdf.sh can then be used to convert the resulting html file into a PDF document (using headless Chrome).
"""

from pathlib import Path

def main():
    # Read the config file to get the list of HTML files to combine
    config_path = Path(__file__).parent / "html2doc.cfg"
    with open(config_path, "r") as f:
        html_files = [line.strip() for line in f if line.strip() and not line.startswith("#")]

    # Read the first HTML file (which will be the base for the combined document)
    combined_html = ""
    for i, html_file in enumerate(html_files):
        html_path = Path(__file__).parent / html_file
        with open(html_path, "r", encoding="utf-8") as f:
            html_content = f.read()
        if i == 0:
            html_content = html_content.split("</body>")[0]

        # Remove headers and footers from all but the first file
        if i > 0:
            # This is a very simple way to remove headers and footers, and may need to be adjusted based on the actual structure of the HTML files
            html_content = html_content.split("<body>")[1].split("</body>")[0]

        combined_html += html_content

    combined_html = combined_html.replace('<a href="index.html">Back to index</a>', '<div style="page-break-after: always;"></div>')
    combined_html += "</body></html>"

    # Write the combined HTML to a new file
    output_path = Path(__file__).parent / "complete_NGSG_doc.html"
    with open(output_path, "w", encoding="utf-8") as f:
        f.write(combined_html)

if __name__ == "__main__":
    main()