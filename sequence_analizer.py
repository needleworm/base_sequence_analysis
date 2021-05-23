"""
  Author    Byunghyun Ban
  Email     halfbottle@sangsang.farm

"""

RNA_dict = {
    "UUU" : "Phe", "UUC" : "Phe", "UUA" : "Leu", "UUG" : "Leu",
    "CUU" : "Leu", "CUC" : "Leu", "CUA" : "Leu", "CUG" : "Leu",
    "AUU" : "Ile", "AUC" : "Ile", "AUA" : "Ile", "AUG" : "Met",
    "GUU" : "Val", "GUC" : "Val", "GUA" : "Val", "GUG" : "Val",
    "UCU" : "Ser", "UCC" : "Ser", "UCA" : "Ser", "UCG" : "Ser",
    "CCU" : "Pro", "CCC" : "Pro", "CCA" : "Pro", "CCG" : "Pro",
    "ACU" : "Thr", "ACC" : "Thr", "ACA" : "Thr", "ACG" : "Thr",
    "GCU" : "Ala", "GCC" : "Ala", "GCA" : "Ala", "GCG" : "Ala",
    "UAU" : "Tyr", "UAC" : "Tyr", "UAA" : "Stop", "UAG" : "Stop",
    "CAU" : "His", "CAC" : "His", "CAA" : "Gln", "CAG" : "Gln",
    "AAU" : "Asn", "AAC" : "Asn", "AAA" : "Lys", "AAG" : "Lys",
    "GAU" : "Asp", "GAC" : "Asp", "GAA" : "Glu", "GAG" : "Glu",
    "UGU" : "Cys", "UGC" : "Cys", "UGA" : "Stop", "UGG" : "Trp",
    "CGU" : "Arg", "CGC" : "Arg", "CGA" : "Arg", "CGG" : "Arg",
    "AGU" : "Ser", "AGC" : "Ser", "AGA" : "Arg", "AGG" : "Arg",
    "GGU" : "Gly", "GGC" : "Gly", "GGA" : "Gly", "GGG" : "Gly"
}

DNA_dict = {
    "TTT" : "Phe", "TTC" : "Phe", "TTA" : "Leu", "TTG" : "Leu",
    "CTT" : "Leu", "CTC" : "Leu", "CTA" : "Leu", "CTG" : "Leu",
    "ATT" : "Ile", "ATC" : "Ile", "ATA" : "Ile", "ATG" : "Met",
    "GTT" : "Val", "GTC" : "Val", "GTA" : "Val", "GTG" : "Val",
    "TCT" : "Ser", "TCC" : "Ser", "TCA" : "Ser", "TCG" : "Ser",
    "CCT" : "Pro", "CCC" : "Pro", "CCA" : "Pro", "CCG" : "Pro",
    "ACT" : "Thr", "ACC" : "Thr", "ACA" : "Thr", "ACG" : "Thr",
    "GCT" : "Ala", "GCC" : "Ala", "GCA" : "Ala", "GCG" : "Ala",
    "TAT" : "Tyr", "TAC" : "Tyr", "TAA" : "Stop", "TAG" : "Stop",
    "CAT" : "His", "CAC" : "His", "CAA" : "Gln", "CAG" : "Gln",
    "AAT" : "Asn", "AAC" : "Asn", "AAA" : "Lys", "AAG" : "Lys",
    "GAT" : "Asp", "GAC" : "Asp", "GAA" : "Glu", "GAG" : "Glu",
    "TGT" : "Cys", "TGC" : "Cys", "TGA" : "Stop", "TGG" : "Trp",
    "CGT" : "Arg", "CGC" : "Arg", "CGA" : "Arg", "CGG" : "Arg",
    "AGT" : "Ser", "AGC" : "Ser", "AGA" : "Arg", "AGG" : "Arg",
    "GGT" : "Gly", "GGC" : "Gly", "GGA" : "Gly", "GGG" : "Gly"
}

Amino_dict = {
    "Ala" : "A", "Asp" : "D", "Arg" : "R", "His" : "H", "Lys" : "K",
    "Glu" : "E", "Ser" : "S", "Thr" : "T", "Asn" : "N", "Gln" : "Q",
    "Cys" : "C", "Gly" : "G", "Pro" : "P", "Val" : "V", "Ile" : "I",
    "Leu" : "L", "Met" : "M", "Phe" : "F", "Tyr" : "Y", "Trp" : "W"
}

Start_codon_DNA = "ATG"
Start_codon_RNA = "AUG"
Stop_codons_DNA = ["TAA", "TAG", "TGA"]
Stop_codons_RNA = ["UAA", "UAG", "UGA"]


class sequence_analizer():
    def __init__(self, sequence_filename, mode="DNA"):
        self.mode = mode
        self.filename = sequence_filename
        self.sequence = ""
        self.make_sequence(sequence_filename)
        self.transcriptable_sequences = []
        self.protein = []

    def make_sequence(self, sequence_filename):
        file = open(sequence_filename)
        for line in file:
            line = line.strip()
            if not line:
                continue
            if ">" in line:
                continue
            self.sequence += line
        file.close()

    def find_all_possible_proteins(self):
        if not self.transcriptable_sequences:
            self.find_genes()
        for seq in self.transcriptable_sequences:
            protein = self.translation(seq)
            if protein not in self.protein:
                self.protein.append(protein)
        return self.protein

    def translation(self, refined_sequence):
        if self.mode == "DNA":
            dict = DNA_dict
        else:
            dict = RNA_dict
        Protein_sequence = []
        Gene_sequence = ""
        while refined_sequence:
            codon = refined_sequence[:3]
            refined_sequence = refined_sequence[3:]
            if codon not in dict:
                return Protein_sequence, Gene_sequence
            amino_acid = dict[codon]
            if amino_acid == "Stop":
                return Protein_sequence, Gene_sequence
            Protein_sequence.append(amino_acid)
            Gene_sequence += codon
        return Protein_sequence, Gene_sequence

    @staticmethod
    def transcription(self, sequence):
        return "U".join(sequence.split("T"))

    @staticmethod
    def reverse_transcription(self, sequence):
        return "T".join(sequence.split("U"))

    def find_genes(self):
        if self.mode == "DNA":
            start_codon = Start_codon_DNA
            stop_codon1, stop_codon2, stop_codon3 = Stop_codons_DNA
        else:
            start_codon = Start_codon_RNA
            stop_codon1, stop_codon2, stop_codon3 = Stop_codons_RNA

        num_start_codons = self.sequence.count(start_codon)
        template_sequence = self.sequence
        for num in range(num_start_codons):
            sequence = template_sequence[template_sequence.index(start_codon):]
            while start_codon in sequence:
                sequence = sequence[sequence.index(start_codon):]
                stop_codon_idx = [sequence.find(stop_codon1), sequence.find(stop_codon2), sequence.find(stop_codon3)]
                for i in range(3):
                    if stop_codon_idx[i]%3:
                        stop_codon_idx[i] = -1

                if sum(stop_codon_idx) == -3:
                    if sequence[:int(len(sequence)/3)*3] not in self.transcriptable_sequences:
                        self.transcriptable_sequences.append(sequence[:int(len(sequence)/3)*3])
                    break

                stop_site = len(self.sequence)
                for idx in stop_codon_idx:
                    if idx == -1:
                        continue
                    stop_site = min(idx, stop_site)
                if sequence[:stop_site] not in self.transcriptable_sequences:
                    self.transcriptable_sequences.append(sequence[:stop_site])
                sequence = sequence[stop_site+3:]
                if len(sequence) < 3:
                    break
            template_sequence = template_sequence[3:]

    def save_result_as_FASTA_format(self):
        self.find_all_possible_transcription()
        new_file_for_gene = open(self.filename.split(".")[0] + "_GENE.fasta", 'w')
        new_file_for_protein = open(self.filename.split(".")[0] + "_PROTEIN.fasta", 'w')

        for i, el in enumerate(self.transcriptable_sequences):
            new_file_for_gene.write(">Sequence " + str(i) + "\n")
            new_file_for_gene.write(self.gene_sequence_into_fasta(el) + "\n\n")

        for i, el in enumerate(self.protein):
            new_file_for_protein.write(">Sequence " + str(i) + "\n")
            new_file_for_protein.write(self.protein_sequence_into_fasta(el) + "\n\n")

        new_file_for_gene.close()
        new_file_for_protein.close()

    def protein_sequence_into_fasta(self, sequence):
        fasta = ""
        for el in sequence:
            fasta += Amino_dict[el]
        retval = ""
        while len(fasta) > 70:
            retval += fasta[:70]
            retval += "\n"
            fasta = fasta[70:]
        retval += fasta
        return retval

    def gene_sequence_into_fasta(self, sequence):
        retval = ""
        while len(sequence) > 70:
            retval += sequence[:70]
            retval += "\n"
            sequence = sequence[70:]
        retval += sequence
        return retval

    def find_all_possible_transcription(self):
        if self.mode == "DNA":
            start_codon = Start_codon_DNA
        else:
            start_codon = Start_codon_RNA

        num_start_codons = self.sequence.count(start_codon)
        template_sequence = self.sequence
        for i in range(num_start_codons):
            template_sequence = template_sequence[template_sequence.index(start_codon):]
            protein, gene = self.translation(template_sequence)
            if gene not in self.transcriptable_sequences:
                self.protein.append(protein)
                self.transcriptable_sequences.append(gene)
            if len(template_sequence) > 3:
                template_sequence = template_sequence[3:]
            else:
                return
