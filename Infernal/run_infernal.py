import argparse
import warnings
import subprocess
import os
import time
from concurrent.futures import ThreadPoolExecutor
from Bio import SeqIO
from shutil import copyfile
import glob


# Ignorar avisos de depreciação
warnings.filterwarnings(action='ignore', category=FutureWarning)
warnings.filterwarnings('ignore')


def remove_equal_sequences(finput, foutput):
    """Remove sequências idênticas de um arquivo fasta e grava o resultado em outro arquivo."""
    dataset = {}
    with open(foutput, 'a') as arq:
        for seq_record in SeqIO.parse(finput, "fasta"):
            name = seq_record.name
            seq = seq_record.seq
            if seq not in dataset.values():
                dataset[name] = seq
                arq.write(f">{name}\n{seq}\n")
            else:
                print(f"Sequência Removida por Ser Idêntica: {name}")


def run_mainv2(file_path, output_path, cpu):
    """Executa mainv2.py como subprocesso."""
    try:
        name = os.path.basename(file_path).split('.')[0]
        directory = os.path.join(output_path, name)

        if not os.path.exists(directory):
            os.makedirs(directory)

        copyfile(file_path, os.path.join(directory, os.path.basename(file_path)))

        env = os.environ.copy()
        env["INFERNAL_NCPU"] = str(cpu)
        env["OMP_NUM_THREADS"] = str(cpu)
        subprocess.run(['cmsearch', '--cut_ga', '--rfam', '--nohmmonly', '--cpu', str(cpu), '--tblout',
                        os.path.join(directory, name + '.tblout'), 'Rfam.cm', os.path.join(directory, os.path.basename(file_path))])
        subprocess.run(['esl-sfetch', '--index', os.path.join(directory, os.path.basename(file_path))])
        subprocess.run(['sh', 'run_infernal.sh', os.path.join(directory, name + '.tblout'),
                        os.path.join(directory, os.path.basename(file_path)), os.path.join(directory, 'seqs.fasta')])
        remove_equal_sequences(os.path.join(directory, 'seqs.fasta'), os.path.join(directory, 'prep_seqs.fasta'))

    except Exception as e:
        print(f"Ocorreu um erro: {e}")


def list_files(folder_path, output_path, max_processes, cpu):
    """Lista arquivos com extensões .fa, .fas ou .fasta na pasta fornecida."""
    try:
        file_extensions = ["*.fa", "*.fas", "*.fasta"]
        files = [file for ext in file_extensions for file in glob.glob(os.path.join(folder_path, ext))]

        if not files:
            print(f"Nenhum arquivo encontrado no caminho: {folder_path}")
            return

        chunks = [files[i:i + max_processes] for i in range(0, len(files), max_processes)]

        for chunk in chunks:
            with ThreadPoolExecutor(max_workers=max_processes) as executor:
                executor.map(run_mainv2, chunk, [output_path]*len(chunk), [cpu]*len(chunk))

    except Exception as e:
        print(f"Ocorreu um erro: {e}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Processa arquivos com extensões .fa, .fas ou .fasta na pasta fornecida.")
    parser.add_argument('-p', '--path', required=True, help='Caminho para a pasta que contém os arquivos .fa, .fas ou .fasta')
    parser.add_argument('-o', '--output', required=True, help='Caminho de saída para a pasta extraída')
    parser.add_argument('--max_processes', type=int, default=12, help='Número máximo de processos para execução paralela.')
    parser.add_argument('--cpu', type=int, default=4, help='Número de CPUs para subprocessos.')
    args = parser.parse_args()

    start_time = time.time()
    list_files(args.path, args.output, args.max_processes, args.cpu)
    end_time = time.time()

    execution_time = end_time - start_time
    print(f"Tempo de execução: {execution_time} segundos")



    #python3 run_infernal.py -p /home/anderson/Documents/PhD\ Files/Genomes-test -o /home/anderson/Documents/PhD\ Files/Genomes-test --max_processes 4 --cpu 3

