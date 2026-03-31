from pypdf import PdfReader, PdfWriter

# ---- Specify PDFs here ----
pdf_files = []

for i in range(0, 100):
    pdf_files.append(f"F_2_5_D2_{i}.pdf")

output_file = "merged_output.pdf"

writer = PdfWriter()

for pdf in pdf_files:
    try:
        reader = PdfReader(pdf)
    except:
        continue
    for page in reader.pages:
        writer.add_page(page)

with open(output_file, "wb") as f:
    writer.write(f)

print(f"Merged {len(pdf_files)} files into '{output_file}'")
