from mysql.connector import connection
import pandas as pd

def query_sql(rna_types):
    conn = connection.MySQLConnection(user='rfamro',
                                      host='mysql-rfam-public.ebi.ac.uk',
                                      database='Rfam', port=4497)

    exclusive_rna_types = set()  # Create a set to store exclusive RNA types

    for rna in rna_types:
        df = pd.read_sql(f"""SELECT f.rfam_acc
                            FROM family f
                            WHERE f.type LIKE '%{rna}%'""", conn)
        
        with open(f'Infernal/accessions/{rna}.txt', 'a') as filename:
            for i in range(len(df)):
                filename.write(str(df.iloc[i, 0]) + '\n')

        for i in range(len(df)):
            exclusive_rna_types.add(df.iloc[i, 0])  # Add RNA types to the set

    conn.close()

    # Query all RNA types
    all_rna_types = set()
    conn = connection.MySQLConnection(user='rfamro',
                                      host='mysql-rfam-public.ebi.ac.uk',
                                      database='Rfam', port=4497)

    df = pd.read_sql(f"""SELECT f.rfam_acc
                        FROM family f""", conn)

    for i in range(len(df)):
        all_rna_types.add(df.iloc[i, 0])

    conn.close()

    # Calculate the exclusive RNA types
    exclusive_rna_types = all_rna_types.difference(exclusive_rna_types)

    # Write exclusive RNA types to a file
    with open('Infernal/accessions/unknown.txt', 'w') as filename:
        for rna in exclusive_rna_types:
            filename.write(rna + '\n')

if __name__ == "__main__":
    rna_types = ['rRNA', 'tRNA', 'sRNA', 'Cis-reg']

    query_sql(rna_types)