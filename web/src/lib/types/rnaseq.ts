export interface GeneStats {
  ensg_id: string;
  gene_symbol: string;
  gene_name: string;
  n_pairs: number;
  mean_log2_control: number;
  mean_log2_treated: number;
  paired_log2_fc: number;
  paired_t_stat: number;
  paired_p_value: number;
  paired_padj_bh: number;
  n_samples: number;
  raw_mean: number;
  raw_var: number;
  log2_mean: number;
  log2_var: number;
  pct_zero: number;
  pct_nonzero: number;
  log2_iqr: number;
  neg_log10_padj: number;
  ma_x_mean_log2: number;
  ma_y_log2_fc: number;
  volcano_x_log2_fc: number;
  volcano_y_neg_log10_padj: number;
  regulation: 'up' | 'down' | 'ns';
  is_significant: boolean;
}