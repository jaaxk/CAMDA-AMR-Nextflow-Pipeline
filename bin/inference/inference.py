import torch
import transformers
import sys
import csv
from Bio import SeqIO
from statistics import mode

def load_model(antibiotic, models_path, device, model_max_length):
    model_path = f"{models_path}/{antibiotic}/best"
    model = transformers.AutoModelForSequenceClassification.from_pretrained(
        model_path,
        num_labels=2,
        trust_remote_code=True,
    )
    model.eval()
    model.to(device)
    
    tokenizer = transformers.AutoTokenizer.from_pretrained(
        model_path,
        model_max_length=model_max_length,
        padding_side="right",
        use_fast=True,
        trust_remote_code=True,
    )

    return model, tokenizer

def run_inference(samples_tsv, model, tokenizer, device, num_hits_norm, model_max_length):

    species_mapping = {
    'klebsiella_pneumoniae': 0,
    'streptococcus_pneumoniae': 1,
    'escherichia_coli': 2,
    'campylobacter_jejuni': 3,
    'salmonella_enterica': 4,
    'neisseria_gonorrhoeae': 5,
    'staphylococcus_aureus': 6,
    'pseudomonas_aeruginosa': 7,
    'acinetobacter_baumannii': 8
    }

    per_accession_preds = []
    batch_size = 32
    with open(samples_tsv, 'r') as f:
        reader = csv.reader(f, delimiter='\t')
        for row in reader:
            accession, genus_species, flanks_path = row
            genus_species = species_mapping[genus_species]
            sequences = []
            num_hits = []
            species = []
            hits = len(sequences)/num_hits_norm
            for record in SeqIO.parse(flanks_path, 'fasta'):
                sequences.append(str(record.seq))
                num_hits.append(hits)
                species.append(genus_species)

            with torch.no_grad():
                per_sequence_preds = []
                for start in range(0, len(sequences), batch_size):
                    batch = sequences[start:start+batch_size]
                    num_hits_batch = num_hits[start:start+batch_size]
                    species_batch = species[start:start+batch_size]
                    enc = tokenizer(batch, padding=True, truncation=True, return_tensors='pt', max_length=model_max_length)

                    inputs = {
                        'input_ids': enc['input_ids'].to(device),
                        'attention_mask': enc['attention_mask'].to(device),
                        'num_hits': torch.tensor(num_hits_batch, dtype=torch.float32).unsqueeze(1).to(device),
                        'species': torch.tensor(species_batch, dtype=torch.long).to(device).unsqueeze(1).to(device),
                    }
                    logits = model(**inputs).logits
                    preds = torch.argmax(logits, dim=-1).cpu().numpy()
                    per_sequence_preds.extend(preds)
                
            per_accession_preds.append((accession, mode(per_sequence_preds)))

    return per_accession_preds
                            
                        
    

def main():
    antibiotic = sys.argv[1]
    samples_tsv = sys.argv[2]
    models_path = sys.argv[3]
    model_max_length = int(sys.argv[4])
    num_hits_norm = sys.argv[5]
    num_hits_norm = float(num_hits_norm)

    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")

    model, tokenizer = load_model(antibiotic, models_path, device, model_max_length)

    results = run_inference(samples_tsv, model, tokenizer, device, num_hits_norm, model_max_length)

    with open(f"{antibiotic}_results.tsv", 'w') as f:
        for accession, pred in results:
            f.write(f"{accession}\t{pred}\n")

if __name__ == '__main__':
    main()
    




    
