# Install Sesame if not already installed
if (!requireNamespace("Sesame", quietly = TRUE)) {
    install.packages("Sesame")
}

# Load Sesame library
library(Sesame)

# Assuming you have a VCF file
vcf_file <- "report/gatk/final.output.vcf"

# Load VCF data
vcf_data <- read_vcf(vcf_file)

# Plot data
sesame_plot <- plot(vcf_data)

# Save or display plot
print(sesame_plot)

