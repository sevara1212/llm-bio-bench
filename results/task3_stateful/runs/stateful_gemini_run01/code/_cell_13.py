# Let's print all 9 columns of cluster_means
pd.set_option('display.max_columns', 15)
pd.set_option('display.width', 1000)
print(cluster_means.astype(float).round(3))