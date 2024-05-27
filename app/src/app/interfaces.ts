export interface PaginationItem {
  title: string;
  page: number;
  class?: string;
}

export interface VirusDiscoveryResult {
  status: string;
  error?: string;
  results?: Array<string>;
}
