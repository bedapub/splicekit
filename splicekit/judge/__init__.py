import os
import sys
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
from scipy.stats.stats import pearsonr
import glob
import math
import splicekit
import plotly.express as px
import pandas
import pybio
import numpy

# this produces a juDGE plot from results_edgeR_genes
def process():

    def plot(comp_name):
        plt.figure()
        sns.set(font_scale=0.7)
        sns.set_style("dark")
        sns.set_style("ticks")
        try:
            fig = sns.scatterplot(data={"gene_logFC":data_x[comp_name], "junction_logFC":data_y[comp_name]}, x="gene_logFC", y="junction_logFC", s=20, alpha=0.7, edgecolor='none', color='#888888')
        except:
            return False
        
        plt.plot([-8, 8], [0, 0], color='#999999', linestyle='--', linewidth=0.3, alpha=0.5)
        plt.plot([0, 0], [-8, 8], color='#999999', linestyle='--', linewidth=0.3, alpha=0.5)

        std_gene = numpy.std(data_x[comp_name])
        std_junction = numpy.std(data_y[comp_name])
        try:
            score = std_gene/std_junction
        except:
            score = float("Inf")

        fig.set(title="{comp_name}, junction logFC vs gene logFC, score={score}".format(score="%.3f" % score, comp_name=comp_name))
        fig.set(xlim=(-8,8))
        fig.set(ylim=(-8,8))
        fig.spines['left'].set_linewidth(0.5)
        fig.spines['left'].set_color('#333333')
        fig.spines['bottom'].set_linewidth(0.5)
        fig.spines['bottom'].set_color('#333333')
        fig.spines['top'].set_linewidth(0.5)
        fig.spines['top'].set_color('#333333')
        fig.spines['right'].set_linewidth(0.5)
        fig.spines['right'].set_color('#333333')
        fig.tick_params(axis='x', colors='#333333', width=0.5)
        fig.tick_params(axis='y', colors='#333333', width=0.5)    
        #sns.despine()

        label_points = set()
        temp = []
        temp_x = data_x[comp_name]
        temp_y = data_y[comp_name]
        temp_hover = data_genes[comp_name]
        for dx, dy, dh in zip(temp_x, temp_y, temp_hover):
            temp.append((dx, dy, dh))
        temp.sort(key = lambda x: x[0])
        for el in temp[:3]:
            label_points.add(el)
        temp.sort(key = lambda x: x[0], reverse=True)
        for el in temp[:3]:
            label_points.add(el)
        temp.sort(key = lambda x: x[1])
        for el in temp[:3]:
            label_points.add(el)
        temp.sort(key = lambda x: x[1], reverse=True)
        for el in temp[:3]:
            label_points.add(el)
        label_points = list(label_points)
        for (dx, dy, dh) in label_points:
            plt.text(dx+0.05, dy+0.05, dh, fontdict=dict(color="red", size=8), bbox=dict(pad=0, facecolor="yellow", alpha=0.5))

        plt.savefig("results/judge/plots/{comp_name}.png".format(comp_name=comp_name), dpi=300)
        plt.close()

        df = pd.DataFrame(list(zip(data_x[comp_name], data_y[comp_name], data_genes[comp_name])), columns =['gene_logFC', 'junction_logFC', 'gene_name'])
        fig = px.scatter(df, x="gene_logFC", y="junction_logFC", hover_data=['gene_name'], title=comp_name)
        fig.update_layout(yaxis_range=[-8,8])
        fig.update_layout(xaxis_range=[-8,8])
        fig.update_layout(title_font_size=11)
        fig.update_layout(font={"size":11})
        fig.update_traces(marker={"size":4})
        fig.update_traces(customdata=[df["gene_name"].to_numpy()], hovertemplate="%{z} (gene)")
        fig.write_html("results/judge/plots/{comp_name}.html".format(comp_name=comp_name), full_html=False, include_plotlyjs="cdn")
        return True

    data_x = {}
    data_y = {}
    data_dpfi = {}
    data_genes = {}

    compound_junction_counts = {}
    present_comps = []
    for (comp_name, _, _, _, _) in splicekit.core.annotation.comparisons:
        print(f"[judge] processing {comp_name}")
        gene_data = {}
        gene_junction_data = {}
        f = open(f"results/results_edgeR_genes_all.tab", "rt")
        header = f.readline().replace("\r", "").replace("\n", "").split("\t")
        r = f.readline()
        while r:
            r = r.replace("\r", "").replace("\n", "").split("\t")
            data = dict(zip(header, r))
            if data["comparison"]!=comp_name:
                r = f.readline()
                continue
            gene_data[data["gene_id"]] = (data["gene_name"], float(data["logFC"]))
            r = f.readline()
        f.close()
        f = open("results/results_edgeR_junctions.tab", "rt")
        header = f.readline().replace("\r", "").replace("\n", "").split("\t")
        r = f.readline()
        while r:
            r = r.replace("\r", "").replace("\n", "").split("\t")
            data = dict(zip(header, r))
            if data["comparison"]!=comp_name:
                r = f.readline()
                continue
            if gene_data.get(data["gene_id"], None)!=None:
                if float(data["fdr"])<splicekit.config.edgeR_FDR_thr:
                    data_dpfi.setdefault(comp_name, []).append(float(data["delta_pfi"]))
                    data_y.setdefault(comp_name, []).append(float(data["logFC"]))
                    data_x.setdefault(comp_name, []).append(gene_data[data["gene_id"]][1])
                    data_genes.setdefault(comp_name, []).append(data["gene_name"])
                    gene_junction_data.setdefault(data["gene_id"], []).append(float(data["logFC"]))
            r = f.readline()
        f.close()

        fdata = open(f"results/judge/data/{comp_name}_data.tab", "wt")
        header = ["gene_id", "gene_name", "gene_logFC", "junctions_logFC"]
        fdata.write("\t".join(header) + "\n")
        for gene_id, (gene_name, gene_fdr) in gene_data.items():
            junction_fdr_list = gene_junction_data.get(gene_name, [])
            compound_junction_counts[comp_name] = compound_junction_counts.setdefault(comp_name, 0) + len(junction_fdr_list)
            row = [gene_id, gene_name, gene_fdr, ",".join([str(el) for el in junction_fdr_list])]
            if len(junction_fdr_list)>0:
                fdata.write("\t".join([str(el) for el in row]) + "\n")
        fdata.close()

        if plot(comp_name):
            present_comps.append(comp_name)

    # score compounds
    f = open("results/judge/scored.tab", "wt")
    f.write("compound_name\tjunctions\tstdev_gene\tstdev_junction\tscore\n")
    results = []
    for comp_name, junctions_count in compound_junction_counts.items():
        if comp_name not in present_comps:
            continue
        std_gene = numpy.std(data_x[comp_name])
        std_junction = numpy.std(data_y[comp_name])
        try:
            score = std_gene/std_junction
        except:
            score = float("Inf")
        results.append((score, junctions_count, comp_name, std_gene, std_junction))
    results.sort()
    sorted_comps = []
    for temp in results:
        sorted_comps.append(temp[2])
        row = [temp[2], temp[1], temp[3], temp[4], temp[0]]
        row = [str(x) for x in row]
        f.write("\t".join(row) + "\n")
    f.close()
    
    make_html(sorted_comps)

# this produces a juDGE plot from results/dge
def process_old():

    def plot(comp_name):
        plt.figure()
        sns.set(font_scale=0.7)
        sns.set_style("dark")
        sns.set_style("ticks")
        try:
            fig = sns.scatterplot(data={"gene_logFC":data_x[comp_name], "junction_logFC":data_y[comp_name]}, x="gene_logFC", y="junction_logFC", s=20, alpha=0.7, edgecolor='none', color='#888888')
        except:
            return False
        
        plt.plot([-8, 8], [0, 0], color='#999999', linestyle='--', linewidth=0.3, alpha=0.5)
        plt.plot([0, 0], [-8, 8], color='#999999', linestyle='--', linewidth=0.3, alpha=0.5)

        std_gene = numpy.std(data_x[comp_name])
        std_junction = numpy.std(data_y[comp_name])
        try:
            score = std_gene/std_junction
        except:
            score = float("Inf")

        fig.set(title="{comp_name}, junction logFC vs gene logFC, score={score}".format(score="%.3f" % score, comp_name=comp_name))
        fig.set(xlim=(-8,8))
        fig.set(ylim=(-8,8))
        fig.spines['left'].set_linewidth(0.5)
        fig.spines['left'].set_color('#333333')
        fig.spines['bottom'].set_linewidth(0.5)
        fig.spines['bottom'].set_color('#333333')
        fig.spines['top'].set_linewidth(0.5)
        fig.spines['top'].set_color('#333333')
        fig.spines['right'].set_linewidth(0.5)
        fig.spines['right'].set_color('#333333')
        fig.tick_params(axis='x', colors='#333333', width=0.5)
        fig.tick_params(axis='y', colors='#333333', width=0.5)    
        #sns.despine()

        label_points = set()
        temp = []
        temp_x = data_x[comp_name]
        temp_y = data_y[comp_name]
        temp_hover = data_genes[comp_name]
        for dx, dy, dh in zip(temp_x, temp_y, temp_hover):
            temp.append((dx, dy, dh))
        temp.sort(key = lambda x: x[0])
        for el in temp[:3]:
            label_points.add(el)
        temp.sort(key = lambda x: x[0], reverse=True)
        for el in temp[:3]:
            label_points.add(el)
        temp.sort(key = lambda x: x[1])
        for el in temp[:3]:
            label_points.add(el)
        temp.sort(key = lambda x: x[1], reverse=True)
        for el in temp[:3]:
            label_points.add(el)
        label_points = list(label_points)
        for (dx, dy, dh) in label_points:
            plt.text(dx+0.05, dy+0.05, dh, fontdict=dict(color="red", size=8), bbox=dict(pad=0, facecolor="yellow", alpha=0.5))

        plt.savefig("results/judge/plots/{comp_name}.png".format(comp_name=comp_name), dpi=300)
        plt.close()

        df = pd.DataFrame(list(zip(data_x[comp_name], data_y[comp_name], data_genes[comp_name])), columns =['gene_logFC', 'junction_logFC', 'gene_name'])
        fig = px.scatter(df, x="gene_logFC", y="junction_logFC", hover_data=['gene_name'], title=comp_name)
        fig.update_layout(yaxis_range=[-8,8])
        fig.update_layout(xaxis_range=[-8,8])
        fig.update_layout(title_font_size=11)
        fig.update_layout(font={"size":11})
        fig.update_traces(marker={"size":4})
        fig.update_traces(customdata=[df["gene_name"].to_numpy()], hovertemplate="%{z} (gene)")
        fig.write_html("results/judge/plots/{comp_name}.html".format(comp_name=comp_name), full_html=False, include_plotlyjs="cdn")
        return True

    data_x = {}
    data_y = {}
    data_dpfi = {}
    data_genes = {}

    compound_junction_counts = {}
    present_comps = []
    for (comp_name, _, _, _, _) in splicekit.core.annotation.comparisons:
        print(f"[judge] processing {comp_name}")
        gene_data = {}
        gene_junction_data = {}
        f = open(f"results/dge/dgeTables/topTable-{comp_name}.txt", "rt")
        header = f.readline().replace("\r", "").replace("\n", "").split("\t")
        r = f.readline()
        while r:
            r = r.replace("\r", "").replace("\n", "").split("\t")
            data = dict(zip(header, r))
            gene_data[data["GeneSymbol"]] = float(data["logFC"])
            r = f.readline()
        f.close()
        f = open("results/results_edgeR_junctions.tab", "rt")
        header = f.readline().replace("\r", "").replace("\n", "").split("\t")
        r = f.readline()
        while r:
            r = r.replace("\r", "").replace("\n", "").split("\t")
            data = dict(zip(header, r))
            if data["comparison"]!=comp_name:
                r = f.readline()
                continue
            if gene_data.get(data["gene_name"], None)!=None:
                if float(data["fdr"])<splicekit.config.edgeR_FDR_thr:
                    data_dpfi.setdefault(comp_name, []).append(float(data["delta_pfi"]))
                    data_y.setdefault(comp_name, []).append(float(data["logFC"]))
                    data_x.setdefault(comp_name, []).append(gene_data[data["gene_name"]])
                    data_genes.setdefault(comp_name, []).append(data["gene_name"])
                    gene_junction_data.setdefault(data["gene_name"], []).append(float(data["logFC"]))
            r = f.readline()
        f.close()

        fdata = open(f"results/judge/data/{comp_name}_data.tab", "wt")
        header = ["gene_name", "gene_logFC", "junctions_logFC"]
        fdata.write("\t".join(header) + "\n")
        for gene_name, gene_fdr in gene_data.items():
            junction_fdr_list = gene_junction_data.get(gene_name, [])
            compound_junction_counts[comp_name] = compound_junction_counts.setdefault(comp_name, 0) + len(junction_fdr_list)
            row = [gene_name, gene_fdr, ",".join([str(el) for el in junction_fdr_list])]
            if len(junction_fdr_list)>0:
                fdata.write("\t".join([str(el) for el in row]) + "\n")
        fdata.close()

        if plot(comp_name):
            present_comps.append(comp_name)

    # score compounds
    f = open("results/judge/scored.tab", "wt")
    f.write("compound_name\tjunctions\tstdev_gene\tstdev_junction\tscore\n")
    results = []
    for comp_name, junctions_count in compound_junction_counts.items():
        if comp_name not in present_comps:
            continue
        std_gene = numpy.std(data_x[comp_name])
        std_junction = numpy.std(data_y[comp_name])
        try:
            score = std_gene/std_junction
        except:
            score = float("Inf")
        results.append((score, junctions_count, comp_name, std_gene, std_junction))
    results.sort()
    sorted_comps = []
    for temp in results:
        sorted_comps.append(temp[2])
        row = [temp[2], temp[1], temp[3], temp[4], temp[0]]
        row = [str(x) for x in row]
        f.write("\t".join(row) + "\n")
    f.close()
    
    make_html(sorted_comps)

def make_html(present_comps):
    html_header = """
<html>
    <title>junctions logFC vs gene logFC</title>
    <link rel="preconnect" href="https://fonts.googleapis.com">
    <link rel="preconnect" href="https://fonts.gstatic.com" crossorigin>
    <link href="https://fonts.googleapis.com/css2?family=Open+Sans:wght@300;400;700&display=swap" rel="stylesheet">
    <link href="fontawesome-free-6.1.2-web/css/all.css" rel="stylesheet">
<body>
"""

    html_code = """
<style>
    * {{
        font-family: "Open Sans", sans-serif;
        font-size: 11px;
        text-color: #f1f1f1;
    }}

    table {{
        border-collapse: inherit;
    }}

    td {{
        font-size: 11px;
        text-align: left;
        font-family: "Open Sans", sans-serif;
        border-bottom: 1px dashed #888888;
        border-right: 1px dashed #cccccc;
        padding-bottom: 20px;
        padding-right: 10px;
    }}

    tr {{
    }}
    
    #tooltip {{
        background-color: #999999;
        color: white;
        padding: 5px 10px;
        border-radius: 4px;
        font-size: 13px;
        display: none;
        -webkit-animation: fadein 0.5s;
    }}

    #tooltip[data-show] {{
        display: block;
    }}

    @-webkit-keyframes fadein {{
        from {{ opacity: 0; }}
        to   {{ opacity: 1; }}
    }}

    #arrow,
    #arrow::before {{
    position: absolute;
    width: 8px;
    height: 8px;
    background: inherit;
    }}

    #arrow {{
    visibility: hidden;
    }}

    #arrow::before {{
    visibility: visible;
    content: '';
    transform: rotate(45deg);
    }}

    #tooltip[data-popper-placement^='top'] > #arrow {{
    bottom: -4px;
    }}

    #tooltip[data-popper-placement^='bottom'] > #arrow {{
    top: -4px;
    }}

    #tooltip[data-popper-placement^='left'] > #arrow {{
    right: -4px;
    }}

    #tooltip[data-popper-placement^='right'] > #arrow {{
    left: -4px;
    }}

</style>
    
<script>
    function copy_text(val) {{
        navigator.clipboard.writeText(val);
    }}

    function show(owner) {{
        const tooltip = document.querySelector('#tooltip');
        const popperInstance = Popper.createPopper(owner, tooltip, {{placement: 'right', modifiers: [{{name: 'offset', options: {{offset: [0, 8],}},}},],}});
        tooltip.setAttribute('data-show', '');
        popperInstance.update();
        setTimeout(function(){{hide()}}, 1000);
    }}

    function hide() {{
        tooltip.removeAttribute('data-show');
        tooltip.setAttribute('opacity', '0');
    }}
</script>

<center>

    <table style="width: 100%">
        {tr_blocks}
    </table>

"""

    html_footer = """
</body>
</html>
    """

    tr_blocks = []
    td_blocks = []
    for comp_name in present_comps:
        td_blocks.append("<td>{comp_name}<br><a href='plots/{image}.png' target=_new><img style='width:300px;' src='plots/{image}.png'></a></td>".format(comp_name=comp_name, image=comp_name))
        if len(td_blocks)==3:
            tr_blocks.append("\n<tr>\n" + "\n".join(td_blocks) + "\n</tr>")
            td_blocks = []

    f = open("results/judge/index.html", "wt")
    f.write(html_header + html_code.format(tr_blocks="\n".join(tr_blocks)) + html_footer)
    f.close()

    f = open("results/judge/index_core.html", "wt")
    f.write(html_code.format(tr_blocks="\n".join(tr_blocks)))
    f.close()

    tr_blocks = []
    td_blocks = []
    for comp_name in present_comps:
        subplot_data = open("results/judge/plots/{comp_name}.html".format(comp_name=comp_name)).read()
        td_blocks.append("<td>{subplot_data}</td>".format(subplot_data=subplot_data))
        if len(td_blocks)==1:
            tr_blocks.append("\n<tr>\n" + "\n".join(td_blocks) + "\n</tr>")
            td_blocks = []

    f = open("results/judge/index_interactive.html", "wt")
    f.write(html_header + html_code.format(tr_blocks="\n".join(tr_blocks)) + html_footer)
    f.close()

    f = open("results/judge/index_interactive_core.html", "wt")
    f.write(html_code.format(tr_blocks="\n".join(tr_blocks)))
    f.close()

def read_results(cutoff = 0.05):
    print("cutoff = ", cutoff)
    results = {}
    f = open("results/results_edgeR_exons.tab", "rt")
    header = f.readline().replace("\r", "").replace("\n", "").split("\t")
    r = f.readline()
    while r:
        r = r.replace("\r", "").replace("\n", "").split("\t")
        data = dict(zip(header, r))
        if float(data["fdr"])<cutoff:
            data_comparison = results.get(data["comparison"], {})
            data_gene = data_comparison.get(data["gene_id"], {})
            data_exons = data_gene.get("exons", [])
            data_exons.append(r)
            data_gene["gene_name"] = data["gene_name"].replace("-", "_").replace(" ", "_")
            data_gene["exons"] = data_exons
            data_comparison[data["gene_id"]] = data_gene
            results[data["comparison"]] = data_comparison
        r = f.readline()
    f.close()

    f = open("results/results_edgeR_junctions.tab", "rt")
    header = f.readline().replace("\r", "").replace("\n", "").split("\t")
    r = f.readline()
    while r:
        r = r.replace("\r", "").replace("\n", "").split("\t")
        data = dict(zip(header, r))
        if float(data["fdr"])<cutoff:
            data_comparison = results.get(data["comparison"], {})
            data_gene = data_comparison.get(data["gene_id"], {})
            data_exons = data_gene.get("junctions", [])
            data_exons.append(r)
            data_gene["gene_name"] = data["gene_name"].replace("-", "_").replace(" ", "_")
            data_gene["junctions"] = data_exons
            data_comparison[data["gene_id"]] = data_gene
            results[data["comparison"]] = data_comparison
        r = f.readline()
    f.close()
    return results

def make_exons(num_exons, up_regulated=[], down_regulated=[]):
    exon_size = num_exons
    exon_spacing = num_exons*3
    data = {}
    data["x"] = []
    data["y"] = []
    data["exons"] = []
    palette = []
    cpos = 1
    for x in range(1, num_exons+1):
        data["x"].append(cpos)
        data["x"].append(cpos + exon_size)
        cpos += exon_size + exon_spacing
        data["y"].append(1)
        data["y"].append(1)
        data["exons"].append(x)
        data["exons"].append(x)
    for x in range(1, num_exons+1):
        if x in up_regulated:
            palette.append("#ff0000")
        elif x in down_regulated:
            palette.append("#0000ff")
        else:
            palette.append("#aaaaaa")
    return data, palette, exon_size, exon_spacing

def make_junctions(junctions, up_regulated=[], down_regulated=[], exon_size=4, exon_spacing=10):
    delta = 0.2
    data = {}
    data["x"] = []
    data["y"] = []
    data["junctions"] = []
    palette = []
    cpos = 1
    for index, (exon_start_index, exon_end_index) in enumerate(junctions):
        junction_index = index+1
        junction_start = exon_start_index*exon_size + (exon_start_index-1)*exon_spacing+1+delta  # first exon end
        junction_end =  (exon_end_index-1)*exon_size + (exon_end_index-1)*exon_spacing+1-delta # second exon start
        junction_middle = (junction_start + junction_end) / 2
        data["x"].append(junction_start)
        data["x"].append(junction_middle)
        data["x"].append(junction_end)
        data["y"].append(1)
        if junction_index in up_regulated:
            data["y"].append(1.05)
        elif junction_index in down_regulated:
            data["y"].append(0.95)
        else:
            data["y"].append(1.05)
        data["y"].append(1)
        data["junctions"].append(junction_index)
        data["junctions"].append(junction_index)
        data["junctions"].append(junction_index)
        if junction_index in up_regulated:
            palette.append("#ff0000")
        elif junction_index in down_regulated:
            palette.append("#0000ff")
        else:
            palette.append("#aaaaaa")
    return data, palette

def splicemap_plot(comparison_name, gene_id, gene_name, data_exons, palette_exons, data_junctions=[], palette_junctions=[]):
    plt.figure()
    sns.set(font_scale=0.7)
    sns.set_style("dark")
    sns.set_style("ticks")
    fig = sns.lineplot(data=data_junctions, x="x", y="y", hue="junctions", color='b', linewidth=2, legend=False, palette=palette_junctions)
    sns.lineplot(data=data_exons, x="x", y="y", hue="exons", color='b', linewidth=5, legend=False, palette=palette_exons)
    fig.set(title=f"{gene_id}, {gene_name}")
    fig.spines['left'].set_linewidth(0.5)
    fig.spines['left'].set_color('#333333')
    fig.spines['bottom'].set_linewidth(0.5)
    fig.spines['bottom'].set_color('#333333')
    fig.tick_params(axis='x', colors='#333333', width=0.5)
    fig.tick_params(axis='y', colors='#333333', width=0.5)    
    fig.set(ylim=(0.8,1.2))
    sns.despine()
    plt.savefig(f"results/splicemap/{comparison_name}/gene_{gene_id}_{gene_name}.png", dpi=300)
    plt.close()
    return f"gene_{gene_id}_{gene_name}.png"

def exon_index_junction(exons, junction_start, junction_stop):
    import itertools
    ids = set()
    for index, (exon_start, exon_stop) in enumerate(exons):
        if exon_start<=junction_start<=exon_stop or exon_start<=junction_stop<=exon_stop:
            ids.add(index+1)
    ids = list(ids)
    result = []
    for L in range(len(ids) + 1):
        for subset in itertools.combinations(ids, L):
            if len(subset)==2:
                subset = list(subset)
                subset.sort()
                result.append(tuple(subset))
    return result

def find_exon(exons, exon_start, exon_stop):
    for index, (temp_start, temp_end) in enumerate(exons):
        if temp_start<=exon_start<=exon_stop<=temp_end:
            return index+1
    return None

def splicemap():

    splicekit.core.features.read_genes_exons()
    cutoff = 0.05
    results = read_results(cutoff=cutoff)

    for comparison in splicekit.core.annotation.comparisons:
        comparison_name = comparison[0]
        if comparison_name not in results:
            continue
        os.system(f"mkdir -p results/splicemap/{comparison_name}")
        f = open(f"results/splicemap/{comparison_name}/index.html", "wt")
        f.write(f"<center><font style='font-size:14px'>FDR threshold for exons and junctions = {cutoff}<br><br>\n")
        sorted_results = []

        for gene_id, data_gene in list(results[comparison_name].items())[:200]:
            gene_name = data_gene["gene_name"]
            sig_exons_ind_up = []
            sig_exons_ind_down = []
            sig_junctions_ind_up = []
            sig_junctions_ind_down = []
            if len(data_gene.get("exons", []))>1:
                print("[splicemap]", comparison_name, splicekit.core.annotation.genes[gene_id]["gene_name"])
                exons = list(splicekit.core.annotation.genes[gene_id]["exons"].keys())
                exons.sort()
                exons = pybio.utils.merge_intervals(exons)
                sig_exons = data_gene.get("exons", [])
                for temp in sig_exons:
                    temp2 = (int(temp[7]), int(temp[8]))
                    logFC = float(temp[-4])
                    if logFC>=0:
                        sig_exons_ind_up.append(find_exon(exons, temp2[0], temp2[1]))
                    else:
                        sig_exons_ind_down.append(find_exon(exons, temp2[0], temp2[1]))

                sig_exons_ind_up = list(set(sig_exons_ind_up))
                sig_exons_ind_down = list(set(sig_exons_ind_down))

                sig_junctions = data_gene.get("junctions", [])
                junctions_up = []
                junctions_down = []
                for temp in sig_junctions:
                    temp2 = (int(temp[7]), int(temp[8]))
                    junctions_exons = exon_index_junction(exons, temp2[0], temp2[1])
                    logFC = float(temp[-4])
                    if logFC>=0:
                        junctions_up += junctions_exons
                    else:
                        junctions_down += junctions_exons

                junctions_up = list(set(junctions_up))
                junctions_down = list(set(junctions_down))

                data_exons, palette_exons, exon_size, exon_spacing = make_exons(len(exons), sig_exons_ind_up, sig_exons_ind_down)
                data_junctions, palette_junctions = make_junctions(junctions_up+junctions_down, range(1, len(junctions_up)+1), range(len(junctions_up)+1, len(junctions_up)+len(junctions_down)+1), exon_size, exon_spacing)
                if len(junctions_up+junctions_down)>0:
                    fname = splicemap_plot(comparison_name, gene_id, gene_name, data_exons, palette_exons, data_junctions, palette_junctions)
                    score = len(junctions_up) + len(junctions_down) + len(sig_exons_ind_up) + len(sig_exons_ind_down)
                    sorted_results.append((score, fname))
            sorted_results.sort(reverse=True)
        count = 0
        for (score, fname) in sorted_results:
            count += 1
            f.write(f"<a href={fname} target=_new><img src={fname} width=400></a>\n")
        f.close()