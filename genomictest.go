package main

import (
	"fmt"
	"math/rand"
	"os"
	"time"
)

func main() {
	rand.Seed(time.Now().UnixNano())

	// Generate sample files
	generateFASTA("sample.fasta")
	generateFASTQ("sample.fastq")
	generateVCF("sample.vcf")
	generateBED("sample.bed")

	fmt.Println("✓ Generated sample.fasta - DNA sequences")
	fmt.Println("✓ Generated sample.fastq - Sequencing reads with quality scores")
	fmt.Println("✓ Generated sample.vcf - Variant calls")
	fmt.Println("✓ Generated sample.bed - Genomic intervals")
}

// Generate FASTA file with sample sequences
func generateFASTA(filename string) {
	f, _ := os.Create(filename)
	defer f.Close()

	sequences := []struct {
		id  string
		seq string
	}{
		{"gene_BRCA1", generateDNASequence(500)},
		{"gene_TP53", generateDNASequence(450)},
		{"gene_EGFR", generateDNASequence(600)},
		{"chromosome_fragment", generateDNASequence(1000)},
		{"mitochondrial_DNA", generateDNASequence(350)},
	}

	for _, s := range sequences {
		f.WriteString(fmt.Sprintf(">%s\n", s.id))
		// Write sequence in 80 character lines
		for i := 0; i < len(s.seq); i += 80 {
			end := i + 80
			if end > len(s.seq) {
				end = len(s.seq)
			}
			f.WriteString(s.seq[i:end] + "\n")
		}
	}
}

// Generate FASTQ file with quality scores
func generateFASTQ(filename string) {
	f, _ := os.Create(filename)
	defer f.Close()

	for i := 1; i <= 10; i++ {
		seq := generateDNASequence(150)
		quality := generateQualityScores(150)

		f.WriteString(fmt.Sprintf("@READ_%d\n", i))
		f.WriteString(seq + "\n")
		f.WriteString("+\n")
		f.WriteString(quality + "\n")
	}
}

// Generate VCF file with variant calls
func generateVCF(filename string) {
	f, _ := os.Create(filename)
	defer f.Close()

	// VCF header
	f.WriteString("##fileformat=VCFv4.2\n")
	f.WriteString("##reference=hg38\n")
	f.WriteString("##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total Depth\">\n")
	f.WriteString("##INFO=<ID=AF,Number=A,Type=Float,Description=\"Allele Frequency\">\n")
	f.WriteString("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n")
	f.WriteString("##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype Quality\">\n")
	f.WriteString("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE1\n")

	// Sample variants
	variants := []struct {
		chr, pos, id, ref, alt string
		qual                   int
		dp                     int
		af                     float64
	}{
		{"chr1", "12345", "rs123456", "A", "G", 99, 45, 0.52},
		{"chr1", "67890", "rs234567", "C", "T", 85, 38, 0.48},
		{"chr2", "23456", "rs345678", "G", "A", 95, 52, 0.55},
		{"chr3", "45678", "rs456789", "T", "C", 78, 42, 0.51},
		{"chr7", "89012", "rs567890", "A", "T", 92, 48, 0.49},
		{"chrX", "34567", "rs678901", "C", "G", 88, 40, 0.53},
	}

	for _, v := range variants {
		f.WriteString(fmt.Sprintf("%s\t%s\t%s\t%s\t%s\t%d\tPASS\tDP=%d;AF=%.2f\tGT:GQ\t0/1:99\n",
			v.chr, v.pos, v.id, v.ref, v.alt, v.qual, v.dp, v.af))
	}
}


func generateBED(filename string) {
	f, _ := os.Create(filename)
	defer f.Close()

	
	f.WriteString("track name=sample_features description=\"Sample genomic features\"\n")

	
	features := []struct {
		chr, start, end, name string
		score                 int
	}{
		{"chr1", "1000", "2000", "exon1", 800},
		{"chr1", "3000", "4500", "exon2", 950},
		{"chr2", "5000", "6200", "promoter", 700},
		{"chr3", "10000", "12000", "enhancer", 850},
		{"chr7", "20000", "21500", "intron", 600},
		{"chrX", "30000", "31000", "UTR", 750},
	}

	for _, feat := range features {
		f.WriteString(fmt.Sprintf("%s\t%s\t%s\t%s\t%d\n",
			feat.chr, feat.start, feat.end, feat.name, feat.score))
	}
}

func generateDNASequence(length int) string {
	bases := []rune{'A', 'C', 'G', 'T'}
	sequence := make([]rune, length)

	for i := 0; i < length; i++ {
		sequence[i] = bases[rand.Intn(4)]
	}

	return string(sequence)
}

func generateQualityScores(length int) string {
	scores := make([]rune, length)

	for i := 0; i < length; i++ {

		quality := rand.Intn(21) + 20

		scores[i] = rune(quality + 33)
	}

	return string(scores)
}
