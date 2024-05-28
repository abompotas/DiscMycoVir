export interface PaginationItem {
  title: string;
  page: number;
  class?: string;
}

export interface VirusDiscoveryResponse {
  status: string;
  error?: string;
  results?: Array<string>;
}
