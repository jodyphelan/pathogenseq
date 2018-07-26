def dict_list2tex(l,columns = None, mappings = {}):
    headings = l[0].keys() if not columns else columns
    rows = []
    header = " & ".join([mappings[x].title() if x in mappings else x.title() for x in headings])+"\\tabularnewline"
    for row in l:
        r = " & ".join(["%.3f" % row[x] if type(row[x])==float else str(row[x]).replace("_", " ") for x in headings])
        rows.append(r)
    column_def = "".join(["l" for _ in headings])
    str_rows = "\\tabularnewline\n".join(rows)+"\\tabularnewline"
    return "\\begin{tabular}{%s}\n%s\n\\hline\n%s\n\\hline\n\end{tabular}" % (column_def,header,str_rows)
