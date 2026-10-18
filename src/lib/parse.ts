// Input parsing + validation for the peptide batch box.

import { AA } from "./features";

export interface ParsedPeptide {
  label: string;
  seq: string;
  error?: string;
}

const ALPHABET = new Set(AA);
const MAX_LEN = 80;

export function parseInput(text: string): ParsedPeptide[] {
  const out: { label: string; seq: string }[] = [];
  let cur: { label: string; seq: string } | null = null;
  const flush = () => {
    if (cur && cur.seq) out.push(cur);
    cur = null;
  };
  for (const rawLine of text.split(/\r?\n/)) {
    const line = rawLine.trim();
    if (!line) continue;
    if (line.startsWith(">")) {
      flush();
      cur = { label: line.slice(1).trim() || `peptide ${out.length + 1}`, seq: "" };
    } else if (cur) {
      cur.seq += line.replace(/[\s,]/g, "");
    } else {
      for (const tok of line.split(/[\s,]+/)) {
        if (tok) out.push({ label: `peptide ${out.length + 1}`, seq: tok });
      }
    }
  }
  flush();
  return out.map((p) => {
    const seq = p.seq.toUpperCase();
    return { ...p, seq, error: validate(seq) };
  });
}

export function validate(seq: string): string | undefined {
  if (seq.length === 0) return "empty sequence";
  if (seq.length > MAX_LEN) return `longer than ${MAX_LEN} residues`;
  const bad = new Set([...seq].filter((c) => !ALPHABET.has(c)));
  if (bad.size) return `non-standard residue(s): ${[...bad].join(", ")}`;
  return undefined;
}
