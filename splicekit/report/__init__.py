import os
import gzip
from datetime import datetime
import splicekit

module_desc = "splicekit | report |"

edgeR_columns_junctions = ["result_id", "comparison", "feature_id", "gene_id", "gene_name", "jbrowse_url", "annotated", "logFC", "fdr", "donor_pattern"]
edgeR_columns_exons = ["result_id", "comparison", "feature_id", "gene_id", "gene_name", "jbrowse_url", "logFC", "fdr"]
edgeR_columns_genes = ["result_id", "comparison", "feature_id", "gene_id", "gene_name", "jbrowse_url", "logFC", "fdr"]

edgeR_results_max = 3000

html_report = """
<html>
    <head>
        <link rel="preconnect" href="https://fonts.googleapis.com">
        <link rel="preconnect" href="https://fonts.gstatic.com" crossorigin>
        <link href="https://fonts.googleapis.com/css2?family=Roboto:wght@100;300;400;500;700&display=swap" rel="stylesheet">        

        <script src="https://code.jquery.com/jquery-3.7.0.js"></script>
        <script src="https://cdn.datatables.net/1.13.7/js/jquery.dataTables.min.js"></script>
        <link rel="stylesheet" type="text/css" href="https://cdn.datatables.net/1.13.7/css/jquery.dataTables.min.css">
    </head>

    <style>
        @import url('https://fonts.googleapis.com/css2?family=Roboto:wght@100;300;400;500;700&display=swap');
    </style>
    
    <style>
        * {{
            font-family: 'Roboto', sans-serif;
            font-size: 12px;
            color: #444444;
        }}

        a {{
            text-decoration: none;
            color: #0000ff;
        }}

        .menu_div {{
            position: fixed;
            top: 2;
            left: 50%;
            color: #444444;
            margin-left: -710px;
            padding-right: 25px;
            padding-left: 5px;
            padding-top: 5px;
            padding-bottom: 5px;
            margin-right: 15px;
            z-index: 1000;
            text-align: left;
            font-size: 13px !important;
            border-right: 1px dashed #cacaca;
            height: 100%;
        }}

        .menu_div a {{
            font-size: 12px !important;
        }}

    </style>

    <body>
    <center>

    <div id="div_splicekit" style="width:1000px; text-align:left">
    <div style='font-size: 13px; color: #8B0000; background-color: #eaeaea; margin-left: -15px; padding-right: 5px; margin-bottom: 15px;'>Project information</div>
    {project_description}
    </div>

    <br>

    <div id="div_A1" style="width:1000px; text-align:left">
    <div style='font-size: 13px; color: #8B0000; background-color: #eaeaea; margin-left: -15px; padding-right: 5px; margin-bottom: 15px;'>Junction level analysis (top {edgeR_results_max})</div>
    Analysis of differential junction usage: create a table of junction read counts for each sample -> then perform edgeR diffSpliceDGE on the count data. For details, check the <a target=_new href='https://github.com/bedapub/splicekit/blob/main/splicekit/core/comps_edgeR.R'>R code of the analysis</a>. 
    <br><br>

    <a href='results/edgeR/junctions_results_fdr005.tab.gz'>Download fdr005 results: junctions_results_fdr005.tab.gz</a>
    <br><br>
    <table id="table_A1" class="display compact" style="width:100%">
        <thead>
        {thead_A1}
        </thead>
        <tbody>
        {tbody_A1}
        </tbody>
        <tfoot>
        {tfoot_A1}
        </tfoot>
    </table>
    </div>

    <br>
    <br>

    <div id="div_A2" style="width:1000px; text-align:left">
    <div style='font-size: 13px; color: #8B0000; background-color: #eaeaea; margin-left: -15px; padding-right: 5px; margin-bottom: 15px;'>Exon level analysis (top {edgeR_results_max})</div>
    Analysis of differential exon usage: create a table of exon read counts for each sample -> then perform edgeR diffSpliceDGE on the count data. For details, check the <a target=_new href='https://github.com/bedapub/splicekit/blob/main/splicekit/core/comps_edgeR.R'>R code of the analysis</a>. 
    <br><br>

    <a href='results/edgeR/exons_results_fdr005.tab.gz'>Download fdr005 results: exons_results_fdr005.tab.gz</a>
    <br><br>
    <table id="table_A2" class="display compact" style="width:100%">
        <thead>
        {thead_A2}
        </thead>
        <tbody>
        {tbody_A2}
        </tbody>
        <tfoot>
        {tfoot_A2}
        </tfoot>
    </table>
    </div>

    <br>
    <br>

    <div id="div_A3" style="width:1000px; text-align:left">
    <div style='font-size: 13px; color: #8B0000; background-color: #eaeaea; margin-left: -15px; padding-right: 5px; margin-bottom: 15px;'>Gene level analysis (top {edgeR_results_max})</div>
    Analysis of differential gene usage: create a table of gene read counts for each sample -> then perform edgeR glmQLFTest on the count data. For details, check the <a target=_new href='https://github.com/bedapub/splicekit/blob/main/splicekit/core/comps_edgeR.R'>R code of the analysis</a>. 
    <br><br>
    
    <a href='results/edgeR/genes_results_fdr005.tab.gz'>Download fdr005 results: genes_results_fdr005.tab.gz</a>
    <br><br>
    <table id="table_A3" class="display compact" style="width:100%">
        <thead>
        {thead_A3}
        </thead>
        <tbody>
        {tbody_A3}
        </tbody>
        <tfoot>
        {tfoot_A3}
        </tfoot>
    </table>
    </div>

    <br><br>

    <div id="div_B" style="width:1000px; text-align:left">
    <div style='font-size: 13px; color: #8B0000; background-color: #eaeaea; margin-left: -15px; padding-right: 5px; margin-bottom: 15px;'>Dispersions</div>
    Dispersion plots (from edgeR analysis) on the level of junction counts, exon counts and gene counts. For details, check the <a target=_new href='https://github.com/bedapub/splicekit/blob/main/splicekit/core/comps_edgeR.R'>R code of the analysis</a>. 
    <br><br>
   
        <center>
        <table border=0>
            <tr>
                <td><a href='results/edgeR/dispersion/junctions_dispersion.png' target=_new><img src='results/edgeR/dispersion/junctions_dispersion.png' width=300></a></td>
                <td><a href='results/edgeR/dispersion/exons_dispersion.png' target=_new><img src='results/edgeR/dispersion/exons_dispersion.png' width=300></a></td>
                <td><a href='results/edgeR/dispersion/genes_dispersion.png' target=_new><img src='results/edgeR/dispersion/genes_dispersion.png' width=300></a></td>
            </tr>
            <tr>
                <td align=center>junctions</td>
                <td align=center>exons</td>
                <td align=center>genes</td>
            </tr>
        </table>
        </center>
    </div>

    <br><br>

    <div id="div_C" style="width:1000px; text-align:left">
    <div style='font-size: 13px; color: #8B0000; background-color: #eaeaea; margin-left: -15px; padding-right: 5px; margin-bottom: 15px;'>juDGE plots</div>
    juDGE (junction differential gene expression) plots. The idea is to plot junction logFC vs. gene logFC for all relevant (FDR<0.05) junctions. A narrow plot (vertical) would suggest splicing changes, vs. a broader plot (horizontal) would suggest more changes on the DGE level.
    <br><br>
    
        <center>
        <table border=0>
            {tbody_C}
        </table>
        </center>
    </div>

    <div class="menu_div">
        <img src='https://raw.githubusercontent.com/bedapub/splicekit/main/media/splicekit_logo.png' width=90>
        <br>
        splicekit {version}<br>
        <a target=_new href='https://github.com/bedapub/splicekit'>GitHub repository</a>
        
        <br><br>

        <a href="#div_splicekit">Project Information</a>
        <br><br>

        <a href="#div_A1">Junctions</a>
        <br>
        <a href="#div_A2">Exons</a>
        <br>
        <a href="#div_A3">Genes</a>
        <br>
        <br>

        <a href="#div_B">Dispersions</a>
        <br>
        <br>

        <a href="#div_C">juDGE plots</a>
        <br>
        <br>
        
        </div>
    
    <script>
        new DataTable('#table_A1');
        new DataTable('#table_A2');
        new DataTable('#table_A3');
    </script>
    
</body>
</html>
"""

def copy_files():
    os.makedirs("report/results/edgeR", exist_ok=True)
    os.makedirs("report/results/edgeR/dispersion", exist_ok=True)
    os.makedirs("report/results/judge/plots", exist_ok=True)

    files = []
    files.append(("results/edgeR/dispersion/*.png", "report/results/edgeR/dispersion"))
    files.append(("results/edgeR/*fdr005*.gz", "report/results/edgeR"))
    files.append(("results/judge/plots/*.html", "report/results/judge/plots"))

    for files_from, files_to in files:
        os.system(f"cp {files_from} {files_to} > /dev/null 2>&1")

def process():
    print(f"{module_desc} generating report")
    copy_files()

    timestamp = datetime.now()
    unique_timestamp_str = timestamp.strftime('%Y%m%d%H%M%S%f')    

    thead_A1 = []
    tbody_A1 = []
    thead_A2 = []
    tbody_A2 = []
    thead_A3 = []
    tbody_A3 = []
    thead_A1.append("<tr>")
    thead_A1.append("\n".join([f"<th>{el}</th>" for el in edgeR_columns_junctions]))
    thead_A1.append("</tr>")
    thead_A2.append("<tr>")
    thead_A2.append("\n".join([f"<th>{el}</th>" for el in edgeR_columns_exons]))
    thead_A2.append("</tr>")
    thead_A3.append("<tr>")
    thead_A3.append("\n".join([f"<th>{el}</th>" for el in edgeR_columns_genes]))
    thead_A3.append("</tr>")

    count = 0
    f = gzip.open("results/edgeR/junctions_results_fdr005.tab.gz", "rt")
    header = f.readline().replace("\r", "").replace("\n", "").split("\t")
    r = f.readline()
    while r:
        r = r.replace("\r", "").replace("\n", "").split("\t")
        data = dict(zip(header, r))
        tbody_A1.append("<tr>")
        data["jbrowse_url"] = f"<a target=_new href={data['jbrowse_url']}>JBrowse2</a>"
        data["logFC"] = format(float(data["logFC"]), ".2e")
        data["fdr"] = format(float(data["fdr"]), ".2e")
        tbody_A1.append("\n".join([f"<td nowrap>{data[el]}</td>" for el in edgeR_columns_junctions]))
        tbody_A1.append("</tr>")
        count += 1
        if count>=edgeR_results_max:
            break
        r = f.readline()
    f.close()

    count = 0
    f = gzip.open("results/edgeR/exons_results_fdr005.tab.gz", "rt")
    header = f.readline().replace("\r", "").replace("\n", "").split("\t")
    r = f.readline()
    while r:
        r = r.replace("\r", "").replace("\n", "").split("\t")
        data = dict(zip(header, r))
        tbody_A2.append("<tr>")
        data["jbrowse_url"] = f"<a target=_new href={data['jbrowse_url']}>JBrowse2</a>"
        data["logFC"] = format(float(data["logFC"]), ".2e")
        data["fdr"] = format(float(data["fdr"]), ".2e")
        tbody_A2.append("\n".join([f"<td nowrap>{data[el]}</td>" for el in edgeR_columns_exons]))
        tbody_A2.append("</tr>")
        count = count + 1
        if count>=edgeR_results_max:
            break
        r = f.readline()
    f.close()

    count = 0
    f = gzip.open("results/edgeR/genes_results_fdr005.tab.gz", "rt")
    header = f.readline().replace("\r", "").replace("\n", "").split("\t")
    r = f.readline()
    while r:
        r = r.replace("\r", "").replace("\n", "").split("\t")
        data = dict(zip(header, r))
        tbody_A3.append("<tr>")
        data["jbrowse_url"] = f"<a target=_new href={data['jbrowse_url']}>JBrowse2</a>"
        data["logFC"] = format(float(data["logFC"]), ".2e")
        data["fdr"] = format(float(data["fdr"]), ".2e")
        tbody_A3.append("\n".join([f"<td nowrap>{data[el]}</td>" for el in edgeR_columns_genes]))
        tbody_A3.append("</tr>")
        count = count + 1
        if count>=edgeR_results_max:
            break
        r = f.readline()
    f.close()

    tbody_C = []
    for comp_name, _, _, _, _ in splicekit.core.annotation.comparisons:
        if os.path.exists(f"results/judge/plots/{comp_name}.html"):
            tbody_C.append(f'<td><object type="text/html" data="results/judge/plots/{comp_name}.html?version={unique_timestamp_str}" width="350" height="400"></object></td>')

    def add_every_n_items(lst, item, n):
        num_insertions = len(lst) // n
        offset = 0
        for i in range(1, num_insertions + 1):
            index = i * n + offset
            lst.insert(index, item)
            offset += 1
        return lst    
    
    tbody_C = add_every_n_items(tbody_C, "</tr><tr>", 3)
    tbody_C = ["<tr>"] + tbody_C
    tbody_C.append("</tr>")

    thead_A1 = "\n".join(thead_A1)
    tbody_A1 = "\n".join(tbody_A1)
    thead_A2 = "\n".join(thead_A2)
    tbody_A2 = "\n".join(tbody_A2)
    thead_A3 = "\n".join(thead_A3)
    tbody_A3 = "\n".join(tbody_A3)
    tbody_C = "\n".join(tbody_C)

    project_descrption = "To display project information, provide a project.description file in the splicekit folder."
    if os.path.exists("project.description"):
        project_descrption = open("project.description").readlines()
        project_descrption = "".join(project_descrption)

    f = open("report/index.html", "wt")
    f.write(html_report.format(version=splicekit.version, edgeR_results_max=edgeR_results_max, project_description=project_descrption, tbody_C=tbody_C, thead_A1=thead_A1, tbody_A1=tbody_A1, tfoot_A1=thead_A1, thead_A2=thead_A2, tbody_A2=tbody_A2, tfoot_A2=thead_A2, thead_A3=thead_A3, tbody_A3=tbody_A3, tfoot_A3=thead_A3))
    f.close()