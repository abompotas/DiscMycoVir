export interface PaginationItem {
  title: string;
  page: number;
  class?: string;
}

export interface VirusDiscoveryResponse {
  status: string;
  error?: string;
  reports?: Array<string>;
  results?: Array<VirusDiscoveryResults>;
}

export interface VirusDiscoveryResults {
  queryName: string;
  queryLetters: number;
  descriptions: Array<VirusDiscoveryResultDescription>;
  alignments: Array<VirusDiscoveryResultAlignment>;
}

export interface VirusDiscoveryResultDescription {
  title: string;
  score: number;
  bits: number;
  e: number;
  numAlignments: number;
}

export interface VirusDiscoveryResultAlignment {
  hitId: string;
  hitDef: string;
  length: number;
  hsps: Array<VirusDiscoveryResultHSP>;
}

export interface VirusDiscoveryResultHSP {
  title?: string;
  length?: number;
  score: number;
  bits: number;
  expect: number;
  numAlignments: number | null;
  identities: number;
  positives: number;
  gaps: number;
  alignLength: number;
  strand: Array<string>;
  frame: Array<number>;
  query: string;
  queryStart: number;
  queryEnd: number;
  match: string;
  sbjct: string;
  sbjctStart: number;
  sbjctEnd: number;
}
