export interface PaginationItem {
  title: string;
  page: number;
  class?: string;
}

export interface Taste {
  umami?: number;
  sweetbitter?: number;
  bitter?: number;
  sweet?: number;
  other?: number;
}

export interface CompoundAnalysis {
  smiles: string;
  structure: string;
  checkad?: boolean;
  taste: Taste;
}

export interface FoodAnalysisResult {
  status: string;
  error?: string;
  result?: Taste;
}

export interface CompoundAnalysisResult {
  status: string;
  error?: string;
  result?: Array<CompoundAnalysis>;
}

export interface FooDBCountResult {
  status: string;
  result?: number;
}

export interface FooDBShort {
  food_id: number
  food_name: string;
  food_name_scientific: string;
  food_group: string;
  food_subgroup: string;
}

export interface FooDBListResult {
  status: string;
  result?: Array<FooDBShort>;
}

export interface FooDBDetailed {
  food_id: number;
  public_id: string;
  food_name: string;
  food_name_scientific: string;
  food_group: string;
  food_subgroup: string;
  food_description: string;
  updated_at: string;
}

export interface FooDBDetailsResult {
  status: string;
  result?: FooDBDetailed;
}

export interface FooDBContentsCountResult {
  status: string;
  result?: number;
}

export interface FooDBContent {
  food_id: number;
  source_id: number;
  name: string | null;
  formula: string | null;
  smiles: string | null;
  orig_content: number | null;
  orig_min: number | null;
  orig_max: number | null;
  orig_unit: string | null;
}

export interface FooDBContentsResult {
  status: string;
  result?: Array<FooDBContent>;
}

export interface FoodTrialResult {
  status: string;
  error?: string;
}

export interface FoodTrialField {
  index: number;
  name: string;
  tasteLabel: string;
  minLabel: string;
  maxLabel: string;
}

export interface PocketomeAnalysisResult {
  status: string;
  results?: PocketomeAnalysis;
}

export interface PocketomeAnalysis {
  results: Array<string>;
  uni_prot_ids: Array<string>;
  graphs: PocketomeAnalysisGraphs;
  centroids: Array<string>;
  chains: PocketomeAnalysisChains;
}

export interface PocketomeAnalysisGraphs {
  filtering: string;
  interactions: string;
  bp_bar: string;
  bp_goterm: string;
  bp_gotree: string;
  cc_bar: string;
  cc_goterm: string;
  cc_gotree: string;
  mf_bar: string;
  mf_goterm: string;
  mf_gotree: string;
  pathway_enrichment: string;
}

export interface PocketomeAnalysisChains {
  protein: string;
  ligand: string;
}
