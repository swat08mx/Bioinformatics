import requests
import numpy as np
from tqdm import tqdm
import pandas as pd
from itertools import combinations
from scipy.stats import mannwhitneyu
from num2words import num2words
from scipy import stats

def fetch_ensembl_id(gene_name):
    url = f"https://rest.ensembl.org/lookup/symbol/homo_sapiens/{gene_name}?content-type=application/json"
    response = requests.get(url)

    if response.status_code == 200:
        data = response.json()
        if 'id' in data:
            return data['id']
        else:
            return None
    else:
        print(f"Failed to retrieve Ensembl ID for {gene_name}. Status code: {response.status_code}")
        return None

#genes = pd.read_csv("/content/drive/MyDrive/final_genes.csv")

selected_genes=[]
# sfari = pd.read_csv("SFARI.csv")
# genes=sfari['gene-symbol'].to_list()
#genes = ["DNAH17","DOCK1","DOCK8","SYNE1","CSMD1","ABCA13","LAMA1","MYO7A","INTS1","NF1","NOTCH1","NRCAM","PCDH15","PCM1","PLCB1","TTN"]
#genes = ["TTN", "WDFY3", "NUP155", "FBN1", "CNOT1", "ABCA13", "MYCBP2", "PRKDC", "HYDIN", "HERC1"]
genes=['TTN']
for gene in tqdm(genes):
  print(f"Currently processing {gene}")
#  gene_name = 'GRIN2A'
  result = fetch_ensembl_id(gene)
  # GRIN2A
  gene_id = result
  url = f"https://rest.ensembl.org/overlap/id/{gene_id}?feature=exon"
  response = requests.get(url, headers={"Content-Type": "application/json"})
  start=[]
  end=[]
  if response.ok:
      exons = response.json()
      for exon in exons:
#          print(f"Exon ID: {exon['id']}, Start: {exon['start']}, End: {exon['end']}, Strand: {exon['strand']}")
          start.append(exon['start'])
          end.append(exon['end'])
  else:
    print(f"Faulty response back")

  dicto={'start':start, 'end':end}
  exon = pd.DataFrame(dicto)
  fixed_log2=[]
  fixed_depth=[]
  names = ["105", "004", "78", "010", "013", "022", "81", "034", "036", "044", "61", "009", "80", "030", "87", "90", "92", "003"]
  for name in names:
    data = pd.read_csv(f"/content/drive/MyDrive/{name}-c.bqsr.cnr", sep="\t")
    data=data.loc[(data['gene']==gene)]
    data.reset_index(drop=False, inplace=True)

    selected_start=[]
    selected_end=[]
    for i in range(len(data)):
      for j in range(len(exon)):
        if exon['start'][j]>= data['start'][i] and exon['end'][j]<= data['end'][i]:
          length_exon = exon['end'][j] - exon['start'][j]
          length_probe = data['end'][i] - data['start'][i]
          overlap_pc = (length_exon/length_probe) * 100
          if overlap_pc >=50.00:
            selected_start.append(data['start'][i])
            selected_end.append(data['end'][i])
#            print(f"Inside overlap is {overlap_pc}")
#         elif exon['start'][j]>= data['start'][i] and exon['end'][j]>= data['end'][i]:
#           length_exon_inside = data['end'][i] - exon['start'][j]
#           length_probe = data['end'][i] - data['start'][i]
#           overlap_pc = (length_exon_inside/length_probe) * 100
# #          print(f"Inside truncated overlap on the right side is {overlap_pc}")
#         elif exon['start'][j]<= data['start'][i] and exon['end'][j]<= data['end'][i]:
#           length_exon_inside = exon['end'][j] - data['start'][i]
#           length_probe = data['end'][i] - data['start'][i]
#           overlap_pc = (length_exon_inside/length_probe) * 100
#          print(f"Inside truncated overlap on the left is {overlap_pc}")
#        elif exon['start'][j]<= data['start'][i] and exon['end'][j]>= data['end'][i]:
#          print(f"Overlap is 100%")
#        else:
#          print("No overlap")
    dicto = {'start':selected_start, 'end':selected_end}
    dats = pd.DataFrame(dicto)
    temp_start=[]
    temp_end=[]
    for i in range(len(dats)):
      if dats['start'][i] not in temp_start:
        temp_start.append(dats['start'][i])
        temp_end.append(dats['end'][i])
    dicto={'start':temp_start, 'end': temp_end}
    dats = pd.DataFrame(dicto)

    data_temp = data[data['start'].isin(dats['start'])]
    data_temp.reset_index(drop=False, inplace=True)

    points=[]
    depth=[]
    for i in range(len(data_temp)):
      #if data['chromosome'][i] == "chr20" and data['start'][i]>=63400208 and data['end'][i]<=63472677:
      points.append(float(f"{data_temp['log2'][i]}"))
      depth.append(data_temp['depth'][i])
    fixed_log2.append(points)
    fixed_depth.append(depth)
    bins=[]
    for i in range(len(points)):
      bins.append(f"Bin{i+1}")
    data = {'Bins':bins, 'log2fc':points}
    df = pd.DataFrame(data)
    import matplotlib.pyplot as plt
    # Plotting the data
    plt.plot(df['Bins'], df['log2fc'], marker='o', linestyle='-', color='b')  # Line plot with dots

  plt.title(f"Log2 Fold Change (log2fc) vs Bins for {gene} gene")
  plt.xlabel("Bins")
  plt.ylabel("Log2 Fold Change (log2fc)")
  plt.grid(True)
  plt.show()
  ## Plot for median curve func
  def median(array, label, tag, gene):
    final_log2 = pd.DataFrame(array)
    final_log2 = final_log2.dropna(axis='columns')
    if tag == "depth":
      plt.figure(figsize=(10, 5))
      plt.boxplot(final_log2)
      plt.title(f"{label} for {gene} gene")
      plt.show()
    temp =[]
    for i in range(len(final_log2.columns)):
      temp.append(num2words(i+1))
    final_log2.columns = temp
    avg = []
    for cols in final_log2.columns:
      avg.append(np.median(np.array(final_log2[cols].to_list())))
    bins = []
    for i in range(len(avg)):
      bins.append(num2words(i+1))
    temp = {'Bins':bins, 'values':avg}
    temp = pd.DataFrame(temp)
    ## code for automatical deletion or amp detection
    pos_counter=0
    neg_counter=0
    for i in range(len(temp['Bins'])):
      if float(temp["values"][i]) > 0.00:
        pos_counter+=1
      else:
        neg_counter+=1
    if pos_counter > neg_counter:
      if tag == "log2":
        plt.figure(figsize=(10, 5))
        plt.xticks(rotation=45)
        plt.plot(temp['Bins'], temp['values'], marker='o', linestyle='-', color='r')
        plt.title(f"{label} for {gene} gene")
    else:
      if tag == "log2":
        plt.figure(figsize=(10, 5))
        plt.xticks(rotation=45)
        plt.plot(temp['Bins'], temp['values'], marker='o', linestyle='-', color='b')
        plt.title(f"{label} for {gene} gene")
    return plt.show()

  result_log2 = median(fixed_log2, "Median curve for the log2fc over the bins", "log2", gene)
  print(result_log2)
  result_depth = median(fixed_depth, "Box plot for depth over the bins", "depth", gene)
  print(result_depth)

  ## code for statistical testing
  final = pd.DataFrame(fixed_log2)
  final = final.dropna(axis='columns')
  index=[]
  for i in range(len(final.columns)):
    index.append(i+1)
  final.columns=index
  import seaborn as sns
  import matplotlib.pyplot as plt
  from matplotlib.colors import LinearSegmentedColormap
  fig, ax = plt.subplots(figsize=(12, 8))
  colors=['blue', 'white', 'red']
  custom_cmap = LinearSegmentedColormap.from_list('custom_cmap', colors)
  sns.heatmap(final, cmap=custom_cmap, center=0.0, robust=True)
  plt.show()
  grp=[]
  groups = list(final.keys())
  for group1, group2 in combinations(groups, 2):
      stat, p = mannwhitneyu(final[group1], final[group2], alternative='two-sided')
      if p <= 0.05:
#        print(f"Comparison {group1} vs {group2}: U={stat}, p={p}, the result is significant in the {gene} gene.")
        grp.append(group1)
        grp.append(group2)
  counts={}
  for i in range(len(grp)):
    if grp[i] not in counts:
      counts[grp[i]]=1
    else:
      counts[grp[i]]=counts[grp[i]]+1
  def bubble_sort_dict_by_values(dictionary):
    # Convert dictionary items to list of tuples
    items_list = list(dictionary.items())

    # Implement bubble sort
    n = len(items_list)
    for i in range(n):
        for j in range(0, n - i - 1):
            if items_list[j][1] > items_list[j + 1][1]:
                items_list[j], items_list[j + 1] = items_list[j + 1], items_list[j]

    # Convert sorted list back to dictionary
    sorted_dict = dict(items_list)

    return sorted_dict

  sorted_dict = bubble_sort_dict_by_values(counts)
  print(sorted_dict)
  names = list(sorted_dict.keys())
  values = list(sorted_dict.values())
  plt.figure(figsize=(12, 6))
  plt.bar(range(len(sorted_dict)), values, tick_label=names)
  plt.show()
  ### standard deviation
  def std_dev(data, tag):
    final = pd.DataFrame(data)
    final = final.dropna(axis='columns')
    index=[]
    for i in range(len(final.columns)):
      index.append(i+1)
    final.columns=index

    sd = np.std(final, axis=0)
    mean = np.mean(final, axis=0)
    cv = ((sd/mean) * 100).to_list()
    nums=[]
    for i in range(len(cv)):
      nums.append(i+1)
    sdf = pd.DataFrame(sd, columns=['std'])
    sdf['index']=nums
    plt.figure(figsize=(10, 5))
    plt.plot(sdf['index'], sdf['std'], marker='o', c='r')
    plt.xticks(ticks=sdf['index'], labels=[str(i) for i in sdf['index']])
    plt.xlabel("Bins")
    plt.ylabel('std')
    plt.title(f"Standard Deviation curve for the {tag} for {gene} gene")
    return plt.show()
  result_log2 = std_dev(fixed_log2, "log2")
  print(result_log2)
  # result_depth = std_dev(fixed_depth, "depth")
  # print(result_depth)

  def analysis(array):
    final_log2 = pd.DataFrame(array)
    final_log2 = final_log2.dropna(axis='columns')
    temp =[]
    for i in range(len(final_log2.columns)):
      temp.append(num2words(i+1))
    final_log2.columns = temp
    avg = []
    for cols in final_log2.columns:
      avg.append(np.median(np.array(final_log2[cols].to_list())))
    bins = []
    for i in range(len(avg)):
      bins.append(num2words(i+1))
    temp = {'Bins':bins, 'values':avg}
    temp = pd.DataFrame(temp)
    sorted = list(map(float, temp['values']))
    sorted.sort()
    return sorted
  res = analysis(fixed_log2)
  print(res)
  ## z score
  z_score = stats.zscore(final, axis=1)
  print("Deletion")
  try:
    for cols in list(z_score.columns):
      dat = z_score.loc[(z_score[cols]<=-3.50)]
      if dat.empty:
        continue
      else:
        print(dat)
        if gene not in selected_genes:
          selected_genes.append(gene)
          print(f"New gene added, The current length of the selection is: {len(selected_genes)}")
  except:
    continue
  print("Amplification")
  try:
    for cols in list(z_score.columns):
      dat = z_score.loc[(z_score[cols]>=3.50)]
      if dat.empty:
        continue
      else:
        print(dat)
        if gene not in selected_genes:
          selected_genes.append(gene)
          print(f"New gene added, The current length of the selection is: {len(selected_genes)}")
  except:
    continue
  ## manual work automated to some extent
  array=fixed_log2
  final_log2 = pd.DataFrame(array)
  final_log2 = final_log2.dropna(axis='columns')
  temp =[]
  sd = np.std(final_log2, axis=0)
  for i in range(len(final_log2.columns)):
    temp.append(num2words(i+1))
  final_log2.columns = temp
  avg = []
  for cols in final_log2.columns:
    avg.append(np.median(np.array(final_log2[cols].to_list())))
  avg = list(map(float, avg))
  bins = []
  for i in range(len(avg)):
    bins.append(num2words(i+1))
  temp = pd.DataFrame(bins, columns=['Bins'])
  temp['values']=avg
  sorted = temp['values'].to_list()
  z_score_log2 = stats.zscore(sorted)
  z_score_log2 = list(map(float, z_score_log2))
  print(z_score_log2)
  z_score_std  = stats.zscore(sd).to_list()
  print(z_score_std)

  bin_pos=[]
  bin_neg=[]
  for i in range(len(z_score_log2)):
    if z_score_log2[i] >= np.quantile(z_score_log2, .50) and z_score_std[i]<=np.quantile(z_score_std, .50):
      bin_pos.append(i+1)
    elif z_score_log2[i] <= np.quantile(z_score_log2, .50) and z_score_std[i]<=np.quantile(z_score_std, .50):
      bin_neg.append(i+1)
  print(f"The positive segment bins are: {bin_pos}")
  print(f"The negative segment bins are: {bin_neg}")
