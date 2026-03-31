from pypdf import PdfReader, PdfWriter

# ---- Specify PDFs here ----
pdf_files = []

for c in ["A", "B"]:
    pdf_files.append(f"{c}1_{c}2.pdf")
    pdf_files.append(f"{c}2_{c}3.pdf")
for i in range(1,4):
    for j in range(1,4):
        pdf_files.append(f"A{i}_B{j}.pdf")

output_file = "merged_output.pdf"

writer = PdfWriter()

for pdf in pdf_files:
    reader = PdfReader(pdf)
    for page in reader.pages:
        writer.add_page(page)

with open(output_file, "wb") as f:
    writer.write(f)

print(f"Merged {len(pdf_files)} files into '{output_file}'")
