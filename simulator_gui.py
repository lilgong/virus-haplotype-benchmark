import tkinter as tk
from tkinter import filedialog, messagebox
import os
import random

transition = {
    'A': 'G',
    'G': 'A',
    'C': 'T',
    'T': 'C',
}

transversion = {
    'A': ['C', 'T'],
    'G': ['C', 'T'],
    'C': ['A', 'G'],
    'T': ['A', 'G'],
}

iupac_codes = {
    "A":["A"],
    "C":["C"],
    "G":["G"],
    "T":["T"],
    "U":["U"],
    "W":["A", "T"],
    "S":["C", "G"],
    "M":["A", "C"],
    "K":["G", "T"],
    "R":["A", "G"],
    "Y":["C", "T"],
    "B":["C", "G", "T"],
    "D":["A", "G", "T"],
    "H":["A", "C", "T"],
    "V":["A", "C", "G"],
    "N":["A", "C", "G", "T"],
    "-":[""]
}


# read fasta file
def read_fasta(file):
    with open(file, 'r') as f:
        lines = f.readlines()
        if not lines or not lines[0].startswith('>'):
            raise ValueError("Input file is not a valid FASTA format: missing header line starting with '>'")
        
        seq = []
        for line in lines[1:]:
            line = line.strip()
            if line.startswith('>'):
                messagebox.showwarning("Warning", "Multiple sequences detected, only use the first sequence.")
                break
            seq.append(line)

    return ''.join(seq)


# validate sequence
def seq_validation(seq):
    seq = seq.upper()
    validated_seq = []
    for base in seq:
        if base in iupac_codes:
            validated_seq.append(random.choice(iupac_codes[base]))
        else:
            raise ValueError(f"Input sequence contains invalid base: {base}")
    return ''.join(validated_seq)


def seq_mutation(seq, mutation_rate, transition_rate=0.7):
    mutated_seq = []
    for base in seq:
        # mutate
        if random.random() <= mutation_rate:
            # transition or transversion
            if random.random() <= transition_rate:
                mutated_seq.append(transition[base])
            else:
                mutated_seq.append(random.choice(transversion[base]))
        # not change
        else:
            mutated_seq.append(base)
    return ''.join(mutated_seq)

def generate_variants(parent_seq, mutation_rate, no_of_variant, mode, transition_rate=0.7):
    variants = []
    # seperate mode
    if mode == "0":
        for i in range(no_of_variant):
            variants.append(seq_mutation(parent_seq, mutation_rate, transition_rate))
    # linked mode
    else:
        child_seq = parent_seq
        for i in range(no_of_variant):
            child_seq = seq_mutation(child_seq, mutation_rate, transition_rate)
            variants.append(child_seq)
    return variants

# build a gui inferface
def gui_interface():
    def select_input():
        input_file = filedialog.askopenfilename(filetypes=[("FASTA files", "*.fasta *.fa"), ("All files", "*.*")])
        if input_file:
            entry_input.delete(0, tk.END)
            entry_input.insert(0, input_file)
            entry_input.xview_moveto(1.0)

    def select_output():
        output_dir = filedialog.askdirectory()
        if output_dir:
            entry_output.delete(0, tk.END)
            entry_output.insert(0, output_dir)
            entry_output.xview_moveto(1.0)

    def start_simulation():
        input_file = entry_input.get().strip()
        output_dir = entry_output.get().strip()
        # if the output_dir does not exist, create one.
        os.makedirs(output_dir, exist_ok=True)
        try:
            mutation_rate = float(entry_mutation.get())
            transition_rate = float(entry_transition.get())
            num_variants = int(entry_num_variants.get())
            mode = int(mode_option.get())
            if not (0 <= mutation_rate <= 1):
                messagebox.showerror("Error", "Mutation rate must be between 0 and 1.")
                return
            if not (0 <= transition_rate <= 1):
                messagebox.showerror("Error", "Transition rate must be between 0 and 1.")
                return
            if num_variants <= 0:
                messagebox.showerror("Error", "Number of variants must be a positive integer.")
                return
        except ValueError:
            messagebox.showerror("Error", "Invalid input")
            return

        if not os.path.isfile(input_file):
            messagebox.showerror("Error", "Cannot find the file.")
            return

        try:
            seq = read_fasta(input_file)
            seq =  seq_validation(seq)
            variants = generate_variants(seq, mutation_rate, num_variants, mode, transition_rate=0.7)
            parent_filename = os.path.join(output_dir, f"parent.fasta")
            with open(parent_filename, "w") as f:
                f.write(f">parent\n{seq}\n")
            for i in range(len(variants)):
                variant_filename = os.path.join(output_dir, f"variant{i+1}.fasta")
                with open(variant_filename,"w") as f:
                    f.write(f">variant{i+1}\n{variants[i]}\n")
            messagebox.showinfo("Success","Simulation done.")
        except Exception as e:
            messagebox.showerror("Error", str(e))

    # build the main window
    window = tk.Tk()
    window.title("Sequence Mutation Simulator")
    window.resizable(False,False)

    # ask user to import files
    tk.Label(window, text="Import your fasta file:").grid(row=0, column=0, sticky="e")
    entry_input = tk.Entry(window, width=30)
    entry_input.grid(row=0, column=1)
    tk.Button(window, text="choose file", command=select_input).grid(row=0, column=2)

    # ask user to give the output file location
    tk.Label(window, text="Output file:").grid(row=1, column=0, sticky="e")
    entry_output = tk.Entry(window, width=30)
    entry_output.insert(0, "variants_output")
    entry_output.grid(row=1, column=1)
    tk.Button(window, text="choose file", command=select_output).grid(row=1, column=2)

    # ask user to specify the mutation rate
    tk.Label(window, text="Mutation rate (0~1):").grid(row=2, column=0, sticky="e")
    entry_mutation = tk.Entry(window, justify="center")
    entry_mutation.insert(0, "0.01")
    entry_mutation.grid(row=2, column=1)

    # ask user to specify the transition rate
    tk.Label(window, text="Transition rate (0~1):").grid(row=3, column=0, sticky="e")
    entry_transition = tk.Entry(window, justify="center")
    entry_transition.insert(0, "0.7")
    entry_transition.grid(row=3, column=1)

    # ask user how many variants need to be generate
    tk.Label(window, text="The number of variants:").grid(row=4, column=0, sticky="e")
    entry_num_variants = tk.Entry(window, justify="center")
    entry_num_variants.insert(0, "5")
    entry_num_variants.grid(row=4, column=1)

    # ask user to select the mode
    tk.Label(window, text="Choose the Mode:").grid(row=5, column=0, sticky="e")
    mode_option = tk.IntVar(value=0)
    tk.Radiobutton(window, text="Separate", variable=mode_option, value=0).grid(row=5, column=1, sticky="w")
    tk.Radiobutton(window, text="Linked", variable=mode_option, value=1).grid(row=5, column=1, sticky="e")

    # start button
    tk.Button(window, text="Start!", command=start_simulation).grid(row=6, column=1, pady=10)

    window.mainloop()

gui_interface()