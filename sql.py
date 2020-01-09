# SQL functionality for SDBL
# Author: Henry Amrhein
# Date: 17 OCT 2019

import gzip
import sqlite3

"""SQL functionality for use by SDBL"""

__version__ = 1.0

TEMP_SCHEMA = "CREATE TEMPORARY TABLE gl(ID INTEGER PRIMARY KEY, name TEXT);"

TEMP_INSERT = "INSERT INTO temp.gl (name) VALUES (?);"

TEMP_DELETE = "DELETE FROM temp.gl;"

ALIAS_SCHEMA = "CREATE TABLE alias (id INTEGER PRIMARY KEY AUTOINCREMENT, prot_id TEXT NOT NULL, alias TEXT NOT NULL, source TEXT NOT NULL);"

ALIAS_INSERT = 'INSERT INTO alias (prot_id, alias, source) VALUES (?, ?, ?);'

ACTIONS_SCHEMA = "CREATE TABLE actions (id INTEGER PRIMARY KEY AUTOINCREMENT, item_id_a TEXT NOT NULL, item_id_b TEXT NOT NULL, mode TEXT NOT NULL, action TEXT, is_directional INT, a_is_acting INT, score INT NOT NULL);"

ACTIONS_INSERT = 'INSERT INTO actions (item_id_a, item_id_b, mode, action, is_directional, a_is_acting, score) VALUES (?, ?, ?, ?, ?, ?, ?);'

EVIDENCE_SCHEMA = "CREATE TABLE evidence (id INTEGER PRIMARY KEY AUTOINCREMENT, protein1 TEXT NOT NULL, protein2 TEXT NOT NULL, neighborhood INT NOT NULL, fusion INT NOT NULL, cooccurence INT NOT NULL, coexpression INT NOT NULL, experimental INT NOT NULL, database INT NOT NULL, textmining INT NOT NULL, combined_score INT NOT NULL);"

EVIDENCE_INSERT = 'INSERT INTO evidence (protein1, protein2, neighborhood, fusion, cooccurence, coexpression, experimental, database, textmining, combined_score) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?);'

ALIAS_QRY = "SELECT prot_id, alias FROM alias JOIN gl ON gl.name = alias.alias;"

ALIAS_REV_QRY = "SELECT prot_id, alias FROM alias JOIN gl ON gl.name = alias.prot_id;"

ACTION_QRY = "SELECT item_id_a, item_id_b, mode, action, is_directional, a_is_acting, score FROM actions INNER JOIN gl ON gl.name = actions.item_id_a AND actions.score > ?;"

EVIDENCE_QRY = 'SELECT protein1, protein2, neighborhood, fusion, cooccurence, coexpression, experimental, database, textmining, combined_score FROM evidence INNER JOIN gl ON gl.name = evidence.protein1 AND combined_score >= ?;'


ACTION_FRAME_COLS = ("gene1", "gene2", "mode", "action", "directional", "gene1_acting", "score")

EVIDENCE_FRAME_COLS = ("neighborhood", "fusion", "cooccurence", "coexpression", "experimental", "database", "textmining", "combined_score")


class SdblSqlException(Exception):
    pass


class SdblSqlCursor:
    def __init__(self, dbh, gene_list):
        self.cursor = dbh.cursor()
        self.gene_set = set(gene_list)

    def __enter__(self):
        self.cursor.executemany(TEMP_INSERT, ((g,) for g in self.gene_set))
        return self.cursor

    def __exit__(self, type, value, traceback):
        self.cursor.execute(TEMP_DELETE)
        self.cursor.close()


class SdblSql:
    def __init__(self, database_file):
        self.dbh = sqlite3.connect(database_file)
        self.cursor = self.dbh.cursor()
        self.cursor.execute(TEMP_SCHEMA)
        self.valid_names = None

    def __del__(self):
        self.dbh.close()

    def close(self):
        self.dbh.close()

    def get_aliases(self, gene_list):
        with SdblSqlCursor(self.dbh, gene_list) as cur:
            cur.execute(ALIAS_QRY)
            aliases = {r[0]: r[1] for r in cur.fetchall()}
        return aliases

    def get_reverse_aliases(self, protein_list, restrict=False):
        with SdblSqlCursor(self.dbh, protein_list) as cur:
            cur.execute(ALIAS_REV_QRY)

            if self.valid_names is not None:
                aliases = {r[0]: r[1] for r in cur.fetchall()
                            if r[1] in self.valid_names}
            else:
                aliases = {r[0]: r[1] for r in cur.fetchall()}

        if restrict:
            aliases = {a for a in aliases}

        return aliases

    def actions_query_gene(self, gene, cutoff_score):
        result = set()
        aliases = self.get_aliases((gene,))

        with SdblSqlCursor(self.dbh, aliases) as cur:
            cur.execute(ACTION_QRY, (cutoff_score,))

            res = [(gene, r[1], r[2], r[3], r[4], r[5], r[6])
                    for r in cur.fetchall()]

        aliases = self.get_reverse_aliases([r[1] for r in res])

        result = {(gene, aliases[r[1]], r[2], r[3], r[4], r[5], r[6]) for
                r in res if r[1] in aliases}

        return sorted(list(result))

    def actions_query_multiple_genes(self, gene_list, cutoff_score):
        result = set()
        aliases = self.get_aliases(gene_list)

        with SdblSqlCursor(self.dbh, aliases) as cur:
            cur.execute(ACTION_QRY, (cutoff_score,))

            res = [(aliases[r[0]], aliases[r[1]], r[2], r[3], r[4], r[5], r[6]) for
                    r in cur.fetchall() if r[1] in aliases]

        return sorted(res)

    def evidence_query_gene(self, gene, cutoff_score):
        aliases = self.get_aliases((gene,))

        with SdblSqlCursor(self.dbh, aliases) as cur:
            cur.execute(EVIDENCE_QRY, (cutoff_score,))

            primary = [[gene] + list(r[1:9]) for r in self.cursor.fetchall()]

        aliases = self.get_reverse_aliases([r[1] for r in primary])

        res = [(r[0], aliases[r[1]], z[0], z[1]) for r in primary
                for z in zip(EVIDENCE_FRAME_COLS, r[2:9]) if z[1]]

        return sorted(res)

    def evidence_query_multiple_genes(self, gene_list, cutoff_score):
        aliases = self.get_aliases(gene_list)

        with SdblSqlCursor(self.dbh, aliases) as cur:
            cur.execute(EVIDENCE_QRY, (cutoff_score,))

            primary = [[aliases[r[0]], aliases[r[1]]] + list(r[2:9]) for r in
                    cur.fetchall() if r[1] in aliases]

        res = [(r[0], r[1], z[0], z[1]) for r in primary
                for z in zip(EVIDENCE_FRAME_COLS, r[2:9]) if z[1]]

        return sorted(res)


def chunked_file(file_handle, rows=10000, delim="\t"):
    """Generator for reading chunks of line-based files."""
    buf = list()
    for line in file_handle:
        if line[0] == 35:
            continue

        if not line:
            break

        tok = line.decode("utf-8").rstrip().split(delim)

        if tok[0] == "protein1" or tok[0] == "item_id_a":
            continue

        if len(tok) == 7:
            tok[4] = 0 if tok[4] == "f" else 1
            tok[5] = 0 if tok[5] == "f" else 1

        buf.append(tok)

        if len(buf) == rows:
            yield buf
            buf = list()

    yield buf

def build_sql_stringdb_database(alias_file, evidence_file, actions_file, database_file, verbose=False):
    """Build Sqlite database from string-db.org organism files"""
    dbh = sqlite3.connect(database_file)
    dbh.isolation_level = None
    cur = dbh.cursor()

    cur.execute("DROP TABLE IF EXISTS alias;")
    cur.execute(ALIAS_SCHEMA)

    if alias_file.endswith("gz"):
        ifs = gzip.open(alias_file, mode="rb")

    else:
        ifs = open(alias_file, mode="r")

    if verbose:
        print("Creating aliases table", end="...", flush=True)

    for chunk in chunked_file(ifs, rows=10000, delim="\t"):

        cur.execute("BEGIN TRANSACTION;")

        cur.executemany(ALIAS_INSERT, chunk)
        cur.execute("COMMIT;")

    ifs.close()

    cur.execute("CREATE INDEX idx_alias ON alias (prot_id, alias);")

    if verbose:
        print("done")

    cur.execute("DROP TABLE IF EXISTS evidence;")
    cur.execute(EVIDENCE_SCHEMA)

    if evidence_file.endswith("gz"):
        ifs = gzip.open(evidence_file, mode="rb")

    else:
        ifs = open(evidence_file, mode="r")

    if verbose:
        print("Creating evidence table", end="...", flush=True)

    for chunk in chunked_file(ifs, rows=10000, delim=" "):
        cur.execute("BEGIN TRANSACTION;")

        cur.executemany(EVIDENCE_INSERT, chunk)
        cur.execute("COMMIT;")

    ifs.close()

    cur.execute("CREATE INDEX idx_evidence ON evidence (protein1, combined_score);")

    if verbose:
        print("done")

    cur.execute("DROP TABLE IF EXISTS actions;")
    cur.execute(ACTIONS_SCHEMA)

    if actions_file.endswith("gz"):
        ifs = gzip.open(actions_file, mode="rb")

    else:
        ifs = open(actions_file, mode="r")

    hdr = ifs.readline()

    if verbose:
        print("Creating molecular action table", end="...", flush=True)

    for chunk in chunked_file(ifs, rows=10000, delim="\t"):
        cur.execute("BEGIN TRANSACTION;")

        cur.executemany(ACTIONS_INSERT, chunk)
        cur.execute("COMMIT;")

    ifs.close()

    cur.execute("CREATE INDEX idx_actions ON actions (item_id_a, score);")

    if verbose:
        print("done")

    dbh.close()

