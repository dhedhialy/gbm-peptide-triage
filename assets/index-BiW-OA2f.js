(function(){const e=document.createElement("link").relList;if(e&&e.supports&&e.supports("modulepreload"))return;for(const s of document.querySelectorAll('link[rel="modulepreload"]'))i(s);new MutationObserver(s=>{for(const n of s)if(n.type==="childList")for(const r of n.addedNodes)r.tagName==="LINK"&&r.rel==="modulepreload"&&i(r)}).observe(document,{childList:!0,subtree:!0});function a(s){const n={};return s.integrity&&(n.integrity=s.integrity),s.referrerPolicy&&(n.referrerPolicy=s.referrerPolicy),s.crossOrigin==="use-credentials"?n.credentials="include":s.crossOrigin==="anonymous"?n.credentials="omit":n.credentials="same-origin",n}function i(s){if(s.ep)return;s.ep=!0;const n=a(s);fetch(s.href,n)}})();const D="ACDEFGHIKLMNPQRSTVWY",le={A:1.8,R:-4.5,N:-3.5,D:-3.5,C:2.5,Q:-3.5,E:-3.5,G:-.4,H:-3.2,I:4.5,L:3.8,K:-3.9,M:1.9,F:2.8,P:-1.6,S:-.8,T:-.7,W:-.9,Y:-1.3,V:4.2},fe={A:71.0788,R:156.1875,N:114.1038,D:115.0886,C:103.1388,Q:128.1307,E:129.1155,G:57.0519,H:137.1411,I:113.1594,L:113.1594,K:128.1741,M:131.1926,F:147.1766,P:97.1167,S:87.0782,T:101.1051,W:186.2132,Y:163.176,V:99.1326},V={D:3.65,E:4.25,C:8.33,Y:10.1},U={K:10.5,R:12,H:6},ye=9.6,ve=2.35,R=7.4,xe=new Set(D.split("").filter(t=>le[t]>0)),$e=["length","net_charge","gravy","aromatic_frac","polar_frac","charged_frac","aliphatic_frac","gly_frac","pro_frac","cys_frac","acidic_frac","basic_frac","max_hydro_run","mw"],z=$e.concat(D.split("").map(t=>`frac_${t}`)),P=(t,e)=>{let a=0;for(const i of t)e.includes(i)&&a++;return a};function q(t){const e=t.length;let a=1/(1+10**(R-ye));a-=1/(1+10**(ve-R));for(const c of Object.keys(V))a-=P(t,c)/(1+10**(V[c]-R));for(const c of Object.keys(U))a+=P(t,c)/(1+10**(R-U[c]));let i=0,s=0;for(const c of t)i=xe.has(c)?i+1:0,s=Math.max(s,i);let n=0,r=0;for(const c of t)n+=le[c],r+=fe[c];const m=[e,a,n/e,P(t,"FWY")/e,P(t,"STNQ")/e,P(t,"DEKRH")/e,P(t,"VLI")/e,P(t,"G")/e,P(t,"P")/e,P(t,"C")/e,P(t,"DE")/e,P(t,"KRH")/e,s,r+18.0153];for(const c of D)m.push(P(t,c)/e);return m}const Be=new Set(D),Q=80;function X(t){const e=[];let a=null;const i=()=>{a&&a.seq&&e.push(a),a=null};for(const s of t.split(/\r?\n/)){const n=s.trim();if(n)if(n.startsWith(">"))i(),a={label:n.slice(1).trim()||`peptide ${e.length+1}`,seq:""};else if(a)a.seq+=n.replace(/[\s,]/g,"");else for(const r of n.split(/[\s,]+/))r&&e.push({label:`peptide ${e.length+1}`,seq:r})}return i(),e.map(s=>{const n=s.seq.toUpperCase();return{...s,seq:n,error:we(n)}})}function we(t){if(t.length===0)return"empty sequence";if(t.length>Q)return`longer than ${Q} residues`;const e=new Set([...t].filter(a=>!Be.has(a)));if(e.size)return`non-standard residue(s): ${[...e].join(", ")}`}function ke(t,e){t.length<e.length&&([t,e]=[e,t]);let a=new Array(e.length+1).fill(0);for(let i=0;i<t.length;i++){const s=[0];for(let n=0;n<e.length;n++)s.push(t[i]===e[n]?a[n]+1:Math.max(a[n+1],s[n]));a=s}return a[e.length]}function Pe(t,e,a=0){if(t===e)return 1;const i=Math.min(t.length,e.length),s=Math.max(t.length,e.length);return 2*i/(i+s)<a?0:2*ke(t,e)/(t.length+e.length)}function ce(t,e,a=3){const i=[];for(let s=0;s<e.length;s++){const n=Pe(t,e[s].s,0);n>0&&i.push({index:s,seq:e[s].s,label:e[s].y,sim:n})}return i.sort((s,n)=>n.sim-s.sim||s.index-n.index),i.slice(0,a)}function _e(t,e){return ce(t,e,1)[0]}function Ae(t,e){let a=0;for(const i of t.trees){let s=0;for(;i.left[s]!==-1;)s=e[i.feature[s]]<=i.thr[s]?i.left[s]:i.right[s];a+=i.prob[s]}return a/t.n_trees}function Ce(t,e){const a=new Map;for(const i of t.trees){let s=0;for(;i.left[s]!==-1;){const n=e[i.feature[s]]<=i.thr[s]?i.left[s]:i.right[s],r=i.feature[s];a.set(r,(a.get(r)??0)+(i.prob[n]-i.prob[s])),s=n}}return[...a.entries()].map(([i,s])=>[i,s/t.n_trees]).filter(([,i])=>Math.abs(i)>1e-9).sort((i,s)=>Math.abs(s[1])-Math.abs(i[1]))}function Me(t,e){return[t<=e+1e-12,1-t<=e+1e-12]}function Se(t){return t[0]&&t[1]?"abstain":t[0]?"negative":t[1]?"positive":"abstain"}function J(t,e,a,i,s,n,r){const m=q(e),c=Ae(t,m),o=Me(c,a),h=ce(e,t.train_seqs,3),p=_e(e,t.train_seqs)?.sim??0,[d,f]=t.train_len_range,b=e.length<d||e.length>f,y=p<i,v=y||b,$=Se(o),k=Ce(t,m).slice(0,4).map(([_,G])=>({feature:s[_],delta:G})),l=`${(100*r).toFixed(0)}%`,g=[];if(y&&g.push(`Outside the training distribution: the closest training peptide is only ${(p*100).toFixed(0)}% similar (flag below the ${(i*100).toFixed(0)}% gate, which sits at the bottom decile of held-out peptide similarity). The ${l} coverage guarantee does not hold here.`),b&&g.push(`Length ${e.length} is outside the ${d}–${f} aa range seen in training, so every descriptor extrapolates. Treated as out-of-domain.`),v)g.push("Verdict withheld: this is an out-of-domain flag, and a conformal guarantee is void outside the distribution it was calibrated on. The probability below is the raw ensemble output, shown for inspection only.");else if($==="abstain"){const _=!o[0]&&!o[1];g.push(_?`Neither label survives the conformal test at the ${l} coverage target: this candidate looks like neither class the calibration set covers. Empty sets are the other edge of the dial — tighten the coverage target far enough and the quantile becomes too strict to admit a label.`:`Both labels survive the conformal test at the ${l} coverage target: the model cannot separate them for this sequence.`)}return!v&&$!=="abstain"&&g.push(`Conformal singleton {${$==="positive"?n:"not "+n}}; ${l} marginal coverage on in-distribution peptides.`),h.length&&g.push(`Closest known peptides: ${h.map(_=>`${_.seq} (${_.label?n:"negative"}, ${(_.sim*100).toFixed(0)}%)`).join("; ")}.`),{prob:c,set:o,q:a,verdict:$,ood:v,maxSim:p,neighbors:h,topAttr:k,reasons:g}}const de=(t,e,a)=>Math.max(e,Math.min(a,t)),Z=new Set("AVILMFWY");function Te(t,e){const a=[];for(let i=0;i+e<=t.length;i++)a.push(t.slice(i,i+e));return a}function Re(t){const e=[],a=Te(t,9);let i;if(a.length>0){let c=0;for(const o of a){const h=[0,3,5,8];c+=h.filter(p=>Z.has(o[p])).length/h.length}i=c/a.length}else i=[...t].filter(o=>Z.has(o)).length/t.length*.5;const s=75*i;e.push({label:"MHC-II anchor density",points:-s,detail:`${(i*100).toFixed(0)}% of 9-mer frames carry a hydrophobic P1/P4/P6/P9 anchor pattern`});const n=t.length>=13&&t.length<=25?1:t.length>=9?.5:.15,r=25*n;e.push({label:"Class-II length window",points:-r,detail:n===1?`${t.length} aa sits in the 13–25 aa optimal MHC-II range`:`${t.length} aa is outside the optimal 13–25 aa range (factor ${n})`});const m=de(100-(s+r),0,100);return{score:Math.round(m),components:e}}const Le=/N[GSTHA]/g,Ee=/D[GST]/g;function Fe(t){const e=[],a=Math.abs(q(t)[1]),i=a>=2&&a<=6?1:a<2?a/2:Math.max(0,1-(a-6)/4);let s=0,n=0;for(const p of t)s="ACILFVM".includes(p)?s+1:0,n=Math.max(n,s);const r=n<=4?1:n>=8?0:1-(n-4)/4,m=100*(.5*i+.5*r);e.push({label:"Solubility (net charge)",points:50*i,detail:`|net charge| = ${a.toFixed(1)} at pH 7.4 (2–6 is the favorable band)`}),e.push({label:"Aggregation (hydrophobic run)",points:50*r,detail:`longest hydrophobic run = ${n} residue${n===1?"":"s"}`});const c=[["Deamidation motifs",S(Le,t),8,"Asn-Gly/Ser/Thr/His/Ala"],["Oxidation (Met)",S(/M/g,t),4,"methionine sulfoxide formation"],["Oxidation (Trp)",S(/W/g,t),6,"tryptophan oxidation"],["Free cysteines",S(/C/g,t),8,"disulfide scrambling / conjugation"],["Asp isomerization motifs",S(Ee,t),6,"Asp-Gly/Ser/Thr"]];let o=0;for(const[p,d,f,b]of c){if(d===0)continue;const y=Math.min(d*f,f*4);o+=y,e.push({label:p,points:-y,detail:`${d} occurrence${d===1?"":"s"} of ${b}`})}t[0]==="Q"&&(o+=6,e.push({label:"N-terminal Gln",points:-6,detail:"pyroglutamate cyclization risk"}));const h=de(m-o,0,100);return{score:Math.round(h),components:e}}function S(t,e){return(e.match(t)??[]).length}const he=[{name:"RGD",pattern:"RGD",target:"integrins αvβ3 / αvβ5",note:"Minimal integrin-recognition motif. Three residues, so it fires often and is not selective on its own; cilengitide is a cyclised RGD pentapeptide that is.",cite:"Nabors et al., Neuro-Oncology 2015 (CORE phase II, cilengitide in GBM)"},{name:"Angiopep-2",pattern:"TFFYGGSRGKRNNFKLKFD",target:"LRP1 (brain transcytosis)",note:"Kunitz-domain peptide that transcytoses across the BBB via LRP1; LRP1 is also detected in glioblastoma biopsies and tumour vasculature.",cite:"Demeule et al., J Neurochem 2008 (doi:10.1111/j.1471-4159.2008.05492.x)"},{name:"TfR-T12",pattern:"THRPPMWSPVWP",target:"transferrin receptor 1",note:"Phage-display TfR1 binder used to target glioma-directed nanoparticles. It binds a receptor that is over-expressed on tumour cells — it is not a demonstrated BBB shuttle.",cite:"Lee et al., Eur J Biochem 2001 (doi:10.1046/j.1432-1327.2001.02073.x)"},{name:"YSA",pattern:"YSAYPDSVPMMS",target:"EphA2",note:"Ephrin-mimetic binder of EphA2, which is up-regulated in tumours and tumour vasculature. A targeting binder, not a BBB-crossing peptide.",cite:"Koolpe, Dail & Pasquale, J Biol Chem 2002 (doi:10.1074/jbc.M208495200)"},{name:"Chlorotoxin",pattern:"MCMPCFTTDHQMARKCDDCCGGKGRGKCYGPQCLCR",target:"MMP-2 / neuropilin-1",note:"36-residue scorpion toxin that binds glioma cells and is the basis of the clinical agent TM-601. Note the 36 aa length: this tool flags it out-of-domain.",cite:"Soroceanu et al., Cancer Res 1998; Deshane et al., J Biol Chem 2003"}],pe=[{target:"EGFRvIII",why:"mutation-specific; addressed with antibodies/vaccines, not a short linear peptide"},{target:"IL13Rα2",why:"glioma-restricted IL-13 receptor; addressed with IL-13 variants and CAR T cells"},{target:"uPAR / PLAUR",why:"invasion and migration marker; addressed with engineered ligands, not plain substrings"},{target:"GD2, B7-H3, HER2",why:"antibody and CAR-T targets in GBM trials; not peptide sequence motifs"}];function Ie(t){const e=[];for(const a of he){let i=0;for(;;){const s=t.indexOf(a.pattern,i);if(s<0)break;e.push({...a,at:s+1}),i=s+1}}return e.sort((a,i)=>a.at-i.at||i.pattern.length-a.pattern.length)}const j=[.4,.3,.2,.15,.1,.05],L=[{key:"bbb",label:"BBB penetration",short:"BBB",source:"trained, conformal"},{key:"cpp",label:"Cell penetration",short:"Cell",source:"trained, conformal"},{key:"immune",label:"Immune compatibility",short:"Immune",source:"heuristic index"},{key:"mfg",label:"Manufacturability",short:"Mfg",source:"heuristic index"}];function N(t,e,a){const i=j[a],s=String(i),n=1-i;return e.map(r=>{if(r.error)return{label:r.label,seq:r.seq,error:r.error};const m=J(t.tasks.bbb,r.seq,t.tasks.bbb.qs[s]??t.tasks.bbb.q,t.ood_sim,t.features,"BBB-penetrating",n),c=J(t.tasks.cpp,r.seq,t.tasks.cpp.qs[s]??t.tasks.cpp.q,t.ood_sim,t.features,"cell-penetrating",n),o=Re(r.seq),h=Fe(r.seq);return{label:r.label,seq:r.seq,bbb:m,cpp:c,immune:o,mfg:h,gbm:Ie(r.seq),objectives:{bbb:m.prob*100,cpp:c.prob*100,immune:o.score,mfg:h.score}}})}const I=`>TfR-T12 (TfR1 binder) — expect confident BBB+
THRPPMWSPVWP
>Penetratin — expect BBB−, cell+ · in the training set
RQIKIWFQNRRMKWKK
>TAT — expect BBB−, cell+
GRKKRRQRRRPPQ
>Angiopep-2 (LRP1 shuttle) — expect abstain
TFFYGGSRGKRNNFKLKFD
>Scrambled Angiopep-2 — expect the same call as its parent
KNKDTKRFYGFGSRLGFFN
>RGD peptide (integrin motif) — expect BBB−
GRGDSPK
>YSA (EphA2 binder) — expect abstain
YSAYPDSVPMMS
>Chlorotoxin (36 aa) — expect out-of-domain
MCMPCFTTDHQMARKCDDCCGGKGRGKCYGPQCLCR
>Held-out random negative — expect confident BBB−
ILACLGFFLDQSVWAH`,ee=[{name:"TfR-T12",expect:"confident BBB+",why:"The only positive call in the batch, and one of the few peptides here that is both novel and in-distribution (closest training peptide 50% similar). This is the model's read, not a literature consensus — the same peptide is a TfR1 binder whose BBB crossing is still debated, which the GBM column below says out loud."},{name:"Penetratin",expect:"BBB−, cell+",why:"Present in the training set, so it is a wiring sanity check rather than a test. A held-out peptide would not come back this clean."},{name:"TAT",expect:"BBB−, cell+",why:"96% similar to a training peptide. Cell-penetrating is not the same as brain-penetrating."},{name:"Angiopep-2",expect:"abstain",why:"The closest training peptide is 79% similar, just under the 80% cluster threshold, and both labels survive the conformal test. Turn the dial to α = 0.30 and it becomes a confident positive call; at α = 0.40 the test admits no label at all, so the set goes empty. This is the coverage/accuracy trade-off, not a bug."},{name:"Scrambled Angiopep-2",expect:"same call and the same probability as its parent",why:"Same residues, shuffled. 33 of the 34 descriptors are permutation-invariant (only the longest hydrophobic run can respond to a shuffle), so the forest cannot see the difference here (measured ΔP = 0.000). The 3-mer baseline in the Validation tab can (ΔP = 0.22) — that is the upgrade path, and the honest ceiling of this feature set."},{name:"RGD peptide",expect:"BBB−, with a GBM lookup hit",why:"Contains the RGD integrin motif. The GBM column is a cited substring lookup, not a binding prediction."},{name:"YSA",expect:"abstain, with a GBM lookup hit",why:"EphA2 binder. No opinion is the honest answer: it is a targeting binder that has never been measured for BBB crossing."},{name:"Chlorotoxin",expect:"out-of-domain",why:"36 aa — outside the 5–30 aa training range, so every descriptor extrapolates and the verdict is withheld. The probability is still shown, raw."},{name:"Held-out random negative",expect:"confident BBB−",why:"A decoy from the negative set that was never trained on, 73% similar to training."}];function u(t){return t.replace(/[&<>"']/g,e=>({"&":"&amp;","<":"&lt;",">":"&gt;",'"':"&quot;","'":"&#39;"})[e])}function te(t,e,a,i){t.querySelectorAll(e).forEach(s=>{s.addEventListener(a,n=>i(s,n))})}function B(t,e=0){return`${(t*100).toFixed(e)}%`}let C="index";function ae(t,e){if(t.ood)return'<span class="chip chip-ood">OOD</span>';const a={positive:["chip chip-pos",`${e}+`],negative:["chip chip-neg",`${e}−`],abstain:["chip chip-abstain","abstain"]},[i,s]=a[t.verdict];return`<span class="${i}" title="${u(t.reasons[0]??"")}">${s}</span>`}function je(t){return`<span class="bar-mini"><i style="width:${(t*100).toFixed(1)}%"></i></span>`}function De(t){return t.length?t.map(e=>`<span class="chip chip-gbm" title="${u(`contains ${e.name} → ${e.target}. ${e.note} ${e.cite}`)}">${u(e.name)}</span>`).join(" "):'<span class="muted">—</span>'}function Ge(t){return`<div class="card inner">
    <h3>GBM-target lookup <span class="tag">cited substring match, not a prediction</span></h3>
    ${t.length?`<table class="mini-table">
        <thead><tr><th>peptide</th><th>target</th><th>at</th><th>what it is</th></tr></thead>
        <tbody>${t.map(a=>`<tr>
              <td class="mono">${u(a.name)}</td>
              <td>${u(a.target)}</td>
              <td class="num">${a.at}</td>
              <td class="hint">${u(a.note)}<br><span class="cite">${u(a.cite)}</span></td>
            </tr>`).join("")}</tbody>
      </table>
      <p class="hint">A hit means the sequence <em>contains</em> a published binder's sequence.
      Affinity, selectivity, stability and BBB crossing are not implied.</p>`:`<p class="hint">No match in this lookup. Targets that have no short linear motif at all —
      ${u(pe.map(a=>a.target).join(", "))} — cannot be covered this way.</p>`}
  </div>`}function Ne(t){if(t.error)return"";const e=q(t.seq),a=n=>e[ue.indexOf(n)]??0,i=t.bbb.topAttr,s=Math.max(...i.map(n=>Math.abs(n.delta)),1e-6);return`
  <tr class="detail-row">
    <td colspan="10">
      <div class="detail-grid">
        <div class="card inner">
          <h3>Why — BBB verdict</h3>
          <p class="prob-line">P(BBB-penetrating) = <strong>${t.bbb.prob.toFixed(3)}</strong>
            &nbsp;·&nbsp; conformal set <code>{${[t.bbb.set[1]?"BBB+":null,t.bbb.set[0]?"BBB−":null].filter(Boolean).join(", ")}}</code>
            &nbsp;·&nbsp; q = ${t.bbb.q.toFixed(3)} &nbsp;·&nbsp; length ${t.seq.length} aa</p>
          <ul class="reasons">${t.bbb.reasons.map(n=>`<li>${u(n)}</li>`).join("")}</ul>
        </div>
        <div class="card inner">
          <h3>Forest attribution (ΔP along the tree paths)</h3>
          <div class="attr">
            ${i.map(n=>{const r=Math.abs(n.delta)/s*50,m=a(n.feature),c=n.delta>=0?50:50-r;return`<div class="attr-row">
                  <span class="attr-name" title="${u(n.feature)}">${u(n.feature)}</span>
                  <span class="attr-track"><i class="attr-fill ${n.delta>=0?"pos":"neg"}" style="left:${c}%;width:${r}%"></i></span>
                  <span class="attr-val">${n.delta>=0?"+":""}${(n.delta*100).toFixed(1)}%</span>
                  <span class="attr-feat">${m.toFixed(2)}</span>
                </div>`}).join("")}
          </div>
          <p class="hint">Bars show how much each descriptor moves the probability
          away from the training prior, averaged over the 200 trees. Right of centre = pushes
          toward BBB+.</p>
        </div>
        <div class="card inner">
          <h3>Closest known peptides <span class="tag">training set</span></h3>
          <table class="mini-table">
            <thead><tr><th>sequence</th><th>label</th><th>similarity</th></tr></thead>
            <tbody>${t.bbb.neighbors.map(n=>`<tr><td class="mono">${u(n.seq)}</td><td>${n.label?"BBB+":"BBB−"}</td><td>${B(n.sim)}</td></tr>`).join("")}</tbody>
          </table>
          <p class="hint">Similarity = 2·LCS/(len<sub>a</sub>+len<sub>b</sub>).
          Below ${B(me)} the peptide is flagged out-of-domain.</p>
        </div>
        <div class="card inner">
          <h3>Cell-penetration verdict</h3>
          <p class="prob-line">P(cell-penetrating) = <strong>${t.cpp.prob.toFixed(3)}</strong></p>
          <ul class="reasons">${t.cpp.reasons.map(n=>`<li>${u(n)}</li>`).join("")}</ul>
        </div>
        ${Ge(t.gbm??[])}
        ${ie("Immune compatibility",t.immune,"risk contribution")}
        ${ie("Manufacturability",t.mfg,"penalty")}
      </div>
    </td>
  </tr>`}function ie(t,e,a){const i=e.components.filter(s=>Math.abs(s.points)>.05).sort((s,n)=>s.points-n.points).map(s=>`<tr>
        <td>${u(s.label)}</td>
        <td class="num ${s.points<0?"neg":"pos"}">${s.points<0?"−":"+"}${Math.abs(s.points).toFixed(0)}</td>
        <td class="hint">${u(s.detail)}${s.points<0?` (${a})`:""}</td>
      </tr>`).join("");return`<div class="card inner">
    <h3>${t} <span class="tag">heuristic index</span></h3>
    <p class="prob-line">score = <strong>${e.score}</strong> / 100</p>
    <table class="mini-table"><tbody>${i}</tbody></table>
  </div>`}let ue=[],me=.4;function Oe(t,e,a){ue=e.model.features,me=e.model.ood_sim;const i=e.results.filter(c=>!c.error),s=i.filter(c=>c.bbb.verdict==="abstain"&&!c.bbb.ood).length,n=i.filter(c=>c.bbb.ood).length,r=i.filter(c=>c.bbb.verdict!=="abstain"&&!c.bbb.ood);t.innerHTML=`
    <section class="card">
      <h2>1 · Paste peptides</h2>
      <textarea id="batch" spellcheck="false" placeholder="One peptide per line, or FASTA:
>name
TFFYGGSRGKRNNFKLKFD">${u(e.input)}</textarea>
      <div class="row">
        <button id="analyze" class="btn btn-primary">Analyze <kbd>⌘⏎</kbd></button>
        <button id="example" class="btn">Load examples</button>
        <button id="clear" class="btn">Clear</button>
      </div>
    </section>
    ${e.input.trim()===I.trim()?qe():""}
    <div class="stats strip">
      ${E("peptides",String(i.length))}
      ${E("confident BBB call",String(r.length))}
      ${E("abstained (BBB)",String(s))}
      ${E("out-of-domain",String(n))}
    </div>

    <section class="card" id="results-card">
      <div class="row row-between">
        <h2>2 · Triage table</h2>
        <label class="inline">sort by
          <select id="sort">
            <option value="index"${C==="index"?" selected":""}>input order</option>
            <option value="bbb"${C==="bbb"?" selected":""}>BBB probability</option>
            <option value="cpp"${C==="cpp"?" selected":""}>cell probability</option>
            <option value="immune"${C==="immune"?" selected":""}>immune compatibility</option>
            <option value="mfg"${C==="mfg"?" selected":""}>manufacturability</option>
          </select>
        </label>
      </div>
      ${e.results.length===0?He():Ke(e)}
    </section>`;const m=t.querySelector("#batch");m.addEventListener("input",()=>a.setInput(m.value)),t.querySelector("#analyze").addEventListener("click",()=>a.analyze()),t.querySelector("#example").addEventListener("click",()=>{a.setInput(I),a.analyze()}),t.querySelector("#clear").addEventListener("click",()=>{a.setInput(""),a.analyze()}),t.querySelector("#sort")?.addEventListener("change",c=>{C=c.target.value,a.reRender()}),te(t,"[data-open]","click",c=>{const o=Number(c.dataset.open);a.select(e.selected===o?null:o)}),te(t,"[data-example]","click",()=>{a.setInput(I),a.analyze()})}function E(t,e){return`<div class="stat"><span class="stat-value">${e}</span><span class="stat-label">${t}</span></div>`}function qe(){return`<details class="card guide" open>
    <summary>What to expect from these ${ee.length} examples</summary>
    <p class="hint">Each peptide shows one behaviour. Read this before the table, so a withheld
    verdict does not read as a failure — abstaining is a designed outcome, not an error.</p>
    <div class="table-wrap"><table class="mini-table">
      <thead><tr><th>peptide</th><th>expected</th><th>why it is in the batch</th></tr></thead>
      <tbody>${ee.map(t=>`<tr>
          <td class="mono">${u(t.name)}</td>
          <td><strong>${u(t.expect)}</strong></td>
          <td class="hint">${u(t.why)}</td>
        </tr>`).join("")}</tbody>
    </table></div>
  </details>`}function He(){return`<div class="empty">
    <p>No batch yet. Paste peptides on the left, or
    <button class="link" data-example>load the example batch</button> —
    literature peptides plus validated B3Pdb positives.</p>
  </div>`}function Ke(t){const e=t.results.map((i,s)=>s);return C!=="index"&&e.sort((i,s)=>{const n=t.results[i].objectives?.[C]??-1;return(t.results[s].objectives?.[C]??-1)-n||i-s}),`<div class="table-wrap"><table class="results">
    <thead><tr>
      <th>#</th><th>name</th><th>sequence</th><th>BBB</th><th>P(BBB)</th>
      <th>cell</th><th title="cited GBM-target lookup — a substring match, not a prediction">GBM</th>
      <th title="0–100, higher is safer (deterministic rule, not a model)">immune<br><small>compat.</small></th>
      <th title="0–100, higher is easier to make (deterministic rule, not a model)">mfg</th>
      <th title="max sequence similarity to any training peptide">sim.</th><th></th>
    </tr></thead>
    <tbody>${e.map(i=>{const s=t.results[i];if(s.error)return`<tr class="row row-error">
          <td class="idx">${i+1}</td>
          <td class="label">${u(s.label)}</td>
          <td class="seq mono">${u(s.seq)}</td>
          <td colspan="8" class="err">skipped — ${u(s.error)}</td>
        </tr>`;const n=t.selected===i;return`<tr class="row${n?" row-open":""}">
        <td class="idx">${i+1}</td>
        <td class="label">${u(s.label)}</td>
        <td class="seq mono">${u(s.seq)}</td>
        <td>${ae(s.bbb,"BBB")}</td>
        <td class="num">${je(s.bbb.prob)} <span class="p">${s.bbb.prob.toFixed(2)}</span></td>
        <td>${ae(s.cpp,"cell")}</td>
        <td>${De(s.gbm??[])}</td>
        <td class="num">${s.immune.score}</td>
        <td class="num">${s.mfg.score}</td>
        <td class="num">${B(s.bbb.maxSim)}</td>
        <td><button class="mini" data-open="${i}" aria-expanded="${n}">${n?"−":"+"}</button></td>
      </tr>
      ${n?Ne(s):""}`}).join("")}</tbody>
  </table></div>`}const Ye="http://www.w3.org/2000/svg";function w(t,e={},a){const i=document.createElementNS(Ye,t);for(const[s,n]of Object.entries(e))i.setAttribute(s,String(n));return a?.appendChild(i),i}function H(t){const e=document.createElement("div");return e.className="chart-tooltip",t.appendChild(e),{show(a,i,s){e.innerHTML=s,e.style.opacity="1",e.style.left=`${a}px`,e.style.top=`${i}px`},hide(){e.style.opacity="0"}}}function K(t,e,a,i){const s=w("svg",{viewBox:`0 0 ${e} ${a}`,class:"chart",style:`max-width:${e}px`});return t.appendChild(s),{svg:s,x:(n,r)=>i.l+(n-r[0])/(r[1]-r[0]||1)*(e-i.l-i.r),y:(n,r)=>a-i.b-(n-r[0])/(r[1]-r[0]||1)*(a-i.t-i.b)}}function Y(t,e,a,i,s,n,r,m,c=o=>String(Math.round(o))){const o=d=>e.l+(d-s[0])/(s[1]-s[0]||1)*(a-e.l-e.r),h=d=>i-e.b-(d-n[0])/(n[1]-n[0]||1)*(i-e.t-e.b);for(let d=0;d<=4;d++){const f=s[0]+(s[1]-s[0])*d/4,b=n[0]+(n[1]-n[0])*d/4;w("line",{x1:e.l,x2:a-e.r,y1:h(b),y2:h(b),class:"grid"},t),w("line",{x1:o(f),x2:o(f),y1:e.t,y2:i-e.b,class:"grid"},t),w("text",{x:o(f),y:i-e.b+16,class:"tick","text-anchor":"middle"},t).textContent=c(f),w("text",{x:e.l-8,y:h(b)+3,class:"tick","text-anchor":"end"},t).textContent=c(b)}w("line",{x1:e.l,x2:a-e.r,y1:i-e.b,y2:i-e.b,class:"axis"},t),w("line",{x1:e.l,x2:e.l,y1:e.t,y2:i-e.b,class:"axis"},t),w("text",{x:(e.l+a-e.r)/2,y:i-6,class:"axis-label","text-anchor":"middle"},t).textContent=r;const p=w("text",{x:12,y:(e.t+i-e.b)/2,class:"axis-label","text-anchor":"middle",transform:`rotate(-90 12 ${(e.t+i-e.b)/2})`},t);p.textContent=m}function We(t,e){t.querySelector("svg")?.remove();const a=e.width??640,i=e.height??420,s={l:46,r:16,t:14,b:40},{svg:n,x:r,y:m}=K(t,a,i,s);Y(n,s,a,i,e.xDomain,e.yDomain,e.xLabel,e.yLabel);const c=H(t);e.points.forEach((h,p)=>{if(h.kind==="ghost")w("circle",{cx:r(h.x,e.xDomain),cy:m(h.y,e.yDomain),r:2.2,class:"pt-ghost"},n);else{const d=w("g",{class:`pt pt-${h.kind}`,tabindex:0},n),f=r(h.x,e.xDomain),b=m(h.y,e.yDomain);w("circle",{cx:f,cy:b,r:h.kind==="front"?5.5:4.5},d),p===e.picked&&w("circle",{cx:f,cy:b,r:9.5,class:"pt-ring"},d);const y=()=>{c.show(f+12,b-8,`<strong>${h.label}</strong>${h.sub?`<br>${h.sub}`:""}`)};d.addEventListener("mouseenter",y),d.addEventListener("focus",y),d.addEventListener("mouseleave",()=>c.hide()),e.onPick&&(d.addEventListener("click",()=>e.onPick(p)),d.addEventListener("keydown",v=>{(v.key==="Enter"||v.key===" ")&&e.onPick(p)}))}});const o=w("rect",{x:s.l,y:s.t,width:a-s.l-s.r,height:i-s.t-s.b,fill:"transparent"},n);o.addEventListener("mousemove",h=>{const p=n.getBoundingClientRect(),d=(h.clientX-p.left)*(a/p.width),f=(h.clientY-p.top)*(i/p.height);let b=-1,y=30;if(e.points.forEach((v,$)=>{if(v.kind==="ghost")return;const k=Math.hypot(r(v.x,e.xDomain)-d,m(v.y,e.yDomain)-f);k<y&&(y=k,b=$)}),b>=0){const v=e.points[b];c.show(r(v.x,e.xDomain)+12,m(v.y,e.yDomain)-8,`<strong>${v.label}</strong>${v.sub?`<br>${v.sub}`:""}`)}else c.hide()}),o.addEventListener("mouseleave",()=>c.hide())}function se(t,e){t.querySelector("svg")?.remove();const a=e.width??420,i=e.height??320,s={l:46,r:16,t:14,b:40},{svg:n,x:r,y:m}=K(t,a,i,s);Y(n,s,a,i,e.xDomain,e.yDomain,e.xLabel,e.yLabel,d=>d.toFixed(2));const c=H(t);e.diagonal&&w("line",{x1:r(e.xDomain[0],e.xDomain),y1:m(e.yDomain[0],e.yDomain),x2:r(e.xDomain[1],e.xDomain),y2:m(e.yDomain[1],e.yDomain),class:"diag"},n),e.series.forEach((d,f)=>{const b=d.points.map((y,v)=>`${v===0?"M":"L"} ${r(y.x,e.xDomain).toFixed(1)} ${m(y.y,e.yDomain).toFixed(1)}`).join(" ");w("path",{d:b,class:`series series-${f}${d.dash?" series-dash":""}`},n),d.markers!==!1&&d.points.forEach((y,v)=>{w("circle",{cx:r(y.x,e.xDomain),cy:m(y.y,e.yDomain),r:v===d.points.length-1?3.8:3.2,class:"marker"},n)})});const o=w("line",{class:"crosshair",y1:s.t,y2:i-s.b,x1:-10,x2:-10},n),h=w("rect",{x:s.l,y:s.t,width:a-s.l-s.r,height:i-s.t-s.b,fill:"transparent"},n);h.addEventListener("mousemove",d=>{const f=n.getBoundingClientRect(),b=(d.clientX-f.left)*(a/f.width),y=e.xDomain[0]+(b-s.l)/(a-s.l-s.r)*(e.xDomain[1]-e.xDomain[0]);o.setAttribute("x1",String(b)),o.setAttribute("x2",String(b));const v=e.series.map($=>{let k=$.points[0];for(const l of $.points)Math.abs(l.x-y)<Math.abs(k.x-y)&&(k=l);return`${$.name}: <strong>${k.y.toFixed(3)}</strong> at ${k.x.toFixed(2)}`}).join("<br>");c.show(b+12,s.t+6,v)}),h.addEventListener("mouseleave",()=>{o.setAttribute("x1","-10"),o.setAttribute("x2","-10"),c.hide()});const p=document.createElement("div");p.className="legend",p.innerHTML=e.series.map(d=>`<span><i class="swatch${d.dash?" swatch-dash":""}"></i>${d.name}</span>`).join(""),t.appendChild(p)}function ne(t,e){t.querySelector("svg")?.remove();const a=e.width??640,i=e.height??320,s={l:46,r:16,t:14,b:62},{svg:n}=K(t,a,i,s);Y(n,s,a,i,[0,e.groups.length],[0,e.max],"",e.yLabel,p=>Number.isInteger(p)?String(p):p.toFixed(1));const r=H(t),m=(a-s.l-s.r)/e.groups.length,c=e.seriesNames.length,o=Math.min(34,(m-14)/c);e.groups.forEach((p,d)=>{const f=s.l+d*m+(m-o*c)/2;p.values.forEach((b,y)=>{const v=(i-s.t-s.b)*b/e.max,$=f+y*o,k=i-s.b-v,l=w("rect",{x:$,y:k,width:o-3,height:v,class:`bar bar-${y}`},n),g=()=>{r.show($+o/2+10,k-6,`<strong>${p.label}</strong><br>${e.seriesNames[y]}: ${b.toFixed(3)}`)};l.addEventListener("mouseenter",g),l.addEventListener("mouseleave",()=>r.hide())}),w("text",{x:s.l+d*m+m/2,y:i-s.b+16,class:"tick","text-anchor":"middle"},n).textContent=p.label,e.annotations?.[d]&&(w("text",{x:s.l+d*m+m/2,y:i-s.b+32,class:"tick tick-note","text-anchor":"middle"},n).textContent=e.annotations[d])});const h=document.createElement("div");h.className="legend",h.innerHTML=e.seriesNames.map((p,d)=>`<span><i class="swatch swatch-${d}"></i>${p}</span>`).join(""),t.appendChild(h)}function oe(t){const e=t.length,a=new Array(e).fill(!0);for(let i=0;i<e;i++)for(let s=0;s<e;s++)i===s||!a[i]||Ve(t[s].dims,t[i].dims)&&(a[i]=!1);return a}function Ve(t,e){let a=!1;for(let i=0;i<t.length;i++){if(t[i]<e[i]-1e-9)return!1;t[i]>e[i]+1e-9&&(a=!0)}return a}let F=null;function Ue(t){if(F)return F;const e=t.model.tasks.bbb.train_seqs.filter(n=>n.y===1),a=Math.max(1,Math.floor(e.length/120)),i=e.filter((n,r)=>r%a===0);return F=N(t.model,i.map((n,r)=>({label:`ref ${r}`,seq:n.s})),4).filter(n=>n.objectives).map(n=>n.objectives),F}function ze(t,e,a){const i=e.results.filter(l=>!l.error&&l.objectives),s=i.map(l=>!l.bbb.ood&&!l.cpp.ood),n=l=>L.find(g=>g.key===l).label,r=l=>L.map(g=>`<option value="${g.key}"${g.key===l?" selected":""}>${g.label} <small>(${g.source})</small></option>`).join(""),m=`
    <div class="row row-between wrap">
      <div class="row wrap">
        <label class="inline">x-axis <select id="px">${r(e.paretoX)}</select></label>
        <label class="inline">y-axis <select id="py">${r(e.paretoY)}</select></label>
        <label class="inline check"><input type="checkbox" id="ghosts"${e.showGhosts?" checked":""}>
          reference cloud <span class="tag">training BBB+</span></label>
      </div>
      <span class="hint">hover to inspect · click a point to pin</span>
    </div>`;if(i.length===0){t.innerHTML=`${m}<div class="empty">
      <p>The Pareto view plots the batch from the <strong>Predict</strong> tab.
      Analyze at least two peptides, then come back — or
      <button class="link" data-goto="predict">go to Predict</button>.</p>
    </div>`,re(t,a);return}const c=l=>[l.objectives.bbb,l.objectives.cpp,l.objectives.immune,l.objectives.mfg],o=s.map((l,g)=>l?g:-1).filter(l=>l>=0),h=oe(o.map(l=>({dims:[i[l].objectives[e.paretoX],i[l].objectives[e.paretoY]]}))),p=i.map(()=>!1);o.forEach((l,g)=>{p[l]=h[g]});const d=oe(o.map(l=>({dims:c(i[l])}))),f=i.map(()=>!1);o.forEach((l,g)=>{f[l]=d[g]});const b=[];if(e.showGhosts)for(const l of Ue(e))b.push({x:l[e.paretoX],y:l[e.paretoY],label:"reference BBB+",sub:"",kind:"ghost"});const y=b.length;i.forEach((l,g)=>{const _=l.objectives,G=p[g];b.push({x:_[e.paretoX],y:_[e.paretoY],label:u(l.label),sub:`${u(l.seq)}<br>${L.map(W=>`${W.short} ${_[W.key].toFixed(0)}`).join(" · ")}${f[g]?"<br><em>on the 4-objective front</em>":""}${s[g]?"":"<br><em>out-of-domain: not ranked</em>"}`,kind:G?"front":!s[g]||l.bbb.verdict==="abstain"?"abstain":"plain"})});const v=i.map((l,g)=>p[g]?{r:l,i:g}:null).filter(Boolean),$=e.picked!==null&&e.picked<i.length?e.picked:null,k=$===null||!s[$]?[]:i.map((l,g)=>({r:l,i:g})).filter(({i:l})=>l!==$&&s[l]&&Qe([i[l].objectives[e.paretoX],i[l].objectives[e.paretoY]],[i[$].objectives[e.paretoX],i[$].objectives[e.paretoY]])).slice(0,4);t.innerHTML=`
    ${m}
    <div class="grid-pareto">
      <div class="card" id="scatter-host"></div>
      <div class="card">
        <h3>Front on the selected axes
          <span class="tag">${v.length} of ${i.length}</span></h3>
        <p class="hint">A peptide is on the front when no other peptide is at least as good on
        both axes and strictly better on one. All four objectives run in the background:
        <strong>${f.filter(Boolean).length} of ${i.length}</strong> peptides are
        non-dominated under all four.${s.includes(!1)?" Out-of-domain peptides (grey, dashed) are plotted but excluded from ranking — their learned objectives are not trustworthy; their heuristic scores still are.":""}</p>
        <div class="chips">
          ${v.map(({r:l,i:g})=>`<button class="chip chip-pick${g===$?" active":""}" data-pick="${g}" title="${u(l.label)}">${u(l.label)}</button>`).join("")||'<span class="hint">no front members</span>'}
        </div>
        ${$===null?'<p class="hint">Pin a point to see its trade-offs here.</p>':`<div class="picked">
                <h3>${u(i[$].label)}</h3>
                <p class="mono">${u(i[$].seq)}</p>
                <table class="mini-table">${L.map(l=>`<tr><td>${l.label}</td><td class="num">${i[$].objectives[l.key].toFixed(0)}</td>
                     <td class="hint">${l.source}</td></tr>`).join("")}</table>
                <p class="hint">${s[$]?k.length===0?"Nothing in the batch beats it on both selected axes.":`On these axes it is dominated by ${k.map(({r:l})=>`<strong>${u(l.label)}</strong>`).join(", ")}.`:"Out-of-domain: shown for inspection, excluded from the front."}</p>
              </div>`}
      </div>
    </div>`,We(t.querySelector("#scatter-host"),{points:b,xLabel:n(e.paretoX),yLabel:n(e.paretoY),xDomain:[0,100],yDomain:[0,100],picked:$===null?void 0:y+$,onPick:l=>{if(l<y)return;const g=l-y;a.pick(e.picked===g?null:g)}}),re(t,a),t.querySelectorAll("[data-pick]").forEach(l=>l.addEventListener("click",()=>a.pick(Number(l.dataset.pick))))}function re(t,e){t.querySelector("#px")?.addEventListener("change",a=>e.setAxis("x",a.target.value)),t.querySelector("#py")?.addEventListener("change",a=>e.setAxis("y",a.target.value)),t.querySelector("#ghosts")?.addEventListener("change",()=>e.toggleGhosts()),t.querySelectorAll("[data-goto]").forEach(a=>a.addEventListener("click",()=>e.goTab("predict")))}function Qe(t,e){let a=!1;for(let i=0;i<t.length;i++){if(t[i]<e[i]-1e-9)return!1;t[i]>e[i]+1e-9&&(a=!0)}return a}function Xe(t,e){const a=e.benchmark;if(!a){t.innerHTML='<div class="empty"><p>Loading benchmark…</p></div>';return}const i=a.b3pred_published_splits.datasets,s={d1:"D1 · B3PP vs CPP",d2:"D2 · vs balanced random",d3:"D3 · vs full random"},n=Object.keys(i).map(o=>({label:s[o]??o,values:[i[o].published_rf_val_auroc,i[o].our_auroc,i[o].our_auroc_novel_only??0]})),r=Object.keys(i).map(o=>`${B(i[o]["frac_valid_ge_0.8_similar_to_train"])} ≥80% like train`),m=[{label:"BBB · random",values:[a.tasks.bbb.random.auroc_mean,a.tasks.bbb.random.auprc_mean]},{label:"BBB · similar",values:[a.tasks.bbb.similarity.auroc_mean,a.tasks.bbb.similarity.auprc_mean]},{label:"Cell · random",values:[a.tasks.cpp.random.auroc_mean,a.tasks.cpp.random.auprc_mean]},{label:"Cell · similar",values:[a.tasks.cpp.similarity.auroc_mean,a.tasks.cpp.similarity.auprc_mean]}],c=(o,h)=>a.tasks[o][h].alpha_sweep.filter(p=>Number.isFinite(p.error_rate)).map(p=>({x:p.abstain_rate,y:p.error_rate}));t.innerHTML=`
    <section class="card">
      <h2>Does it generalise? — the split matters more than the model</h2>
      <p class="lede">B3Pred is the closest published predictor on this exact data. We trained
      our model on <em>their</em> published training partitions and tested on <em>their</em>
      validation partitions, so the test peptides are identical. Then we measured how many of
      those validation peptides are already near-duplicates of training peptides
      (≥80% sequence similarity).</p>
      <div id="bars-a"></div>
      <p class="hint">Note how D3 (full random decoys, 10:1 imbalance) keeps AUROC 0.92 but
      its AUPRC drops to ${i.d3.our_auprc.toFixed(2)} — AUROC hides class imbalance, which is
      exactly why the abstention accounting above is the more useful number.</p>
      <p class="hint">${u(a.b3pred_published_splits.citation)} — published RF AUROC is
      quoted from their paper (their random 80:20 split); “ours” and “ours, novel-only” are
      measured here. Novel-only = validation peptides below 80% similarity to every training
      peptide (n = ${Object.keys(i).map(o=>i[o].n_novel_only).join(" / ")}).</p>
      <div class="table-wrap">
      <table class="mini-table wide">
        <thead><tr><th>split</th><th>published RF</th><th>ours (same split)</th>
        <th>ours, novel-only</th><th>valid ≥80% like train</th><th>valid ≥60% like train</th></tr></thead>
        <tbody>${Object.keys(i).map(o=>`<tr><td>${u(s[o]??o)}</td>
          <td class="num">${i[o].published_rf_val_auroc.toFixed(2)}</td>
          <td class="num">${i[o].our_auroc.toFixed(3)}</td>
          <td class="num">${i[o].our_auroc_novel_only?.toFixed(3)??"—"}</td>
          <td class="num">${B(i[o]["frac_valid_ge_0.8_similar_to_train"],1)}</td>
          <td class="num">${B(i[o]["frac_valid_ge_0.6_similar_to_train"],1)}</td></tr>`).join("")}</tbody>
      </table>
      </div>
    </section>

    <section class="card">
      <h2>Our own protocol: random split vs 80%-similarity cluster split</h2>
      <p class="lede">Random 80:20 splits leak near-duplicates between train and test, as the
      audit above shows. Clustering every peptide at 80% similarity and holding out whole
      clusters gives an honest estimate. Mean ± sd over ${a.tasks.bbb.random.n_seeds} seeds.</p>
      <div id="bars-b"></div>
      <p class="hint">BBB similarity split: ${a.tasks.bbb.similarity.n_clusters} clusters;
      test-set similarity to train pool — median ${B(a.tasks.bbb.similarity.test_to_train_sim.median)},
      ${B(a.tasks.bbb.similarity.test_to_train_sim.frac_below_ood,1)} of test peptides below the
      ${B(a.similarity.ood_threshold)} out-of-domain flag.</p>
    </section>

    <section class="card">
      <h2>The abstention dial, measured</h2>
      <p class="lede">Every point is a different coverage target α, computed on the held-out
      test clusters. Abstaining is not failure — it is the knob that buys a bounded error rate
      on the answers that remain.</p>
      <div class="grid-2">
        <div id="line-abstain"></div>
        <div id="line-reliability"></div>
      </div>
      <p class="hint">Left: error rate among <em>answered</em> peptides vs how many are left
      unanswered. The curve is not monotone step to step — with ~150 calibration peptides each
      α lands on a different empirical quantile, and the middle of the grid is noisy. The ends
      are the reliable statement: refusing ~60% of peptides instead of ~30% roughly halves the
      error rate. Right: reliability diagram (predicted probability vs observed frequency,
      diagonal = perfect calibration); the ~5-bin histogram is coarse for the same
      small-sample reason.</p>
    </section>

    <section class="card">
      <h2>Conformal accounting</h2>
      <div class="table-wrap">
      <table class="mini-table wide">
        <thead><tr><th>task</th><th>split</th><th>coverage (target ${B(1-a.alpha)})</th>
        <th>90% CI</th><th>abstained</th><th>accuracy when answering</th><th>test n</th></tr></thead>
        <tbody>
        ${["bbb","cpp"].flatMap(o=>["random","similarity"].map(h=>{const p=a.tasks[o][h].conformal,d=a.tasks[o][h].coverage_ci;return`<tr><td>${o==="bbb"?"BBB penetration":"Cell penetration"}</td>
              <td>${h}</td><td class="num">${B(p.coverage,1)}</td>
              <td class="num">${d?`[${B(d[0],1)}, ${B(d[1],1)}]`:"—"}</td>
              <td class="num">${B(p.abstain_rate,1)}</td>
              <td class="num">${Number.isFinite(p.singleton_acc)?B(p.singleton_acc,1):"—"}</td>
              <td class="num">${a.tasks[o][h].n_test}</td></tr>`})).join("")}
        </tbody>
      </table>
      </div>
      <p class="hint">Coverage is empirical on the test clusters. The CI is a Wilson interval on
      that empirical coverage — it answers "is the shortfall from the ${B(1-a.alpha)} target
      real, or small-sample noise?", and on ~150 test peptides it is wide by construction.
      The split-conformal guarantee holds marginally on exchangeable data (Lei et al., 2018), and
      is explicitly voided for peptides flagged out-of-domain.</p>
    </section>

    <section class="card">
      <h2>Is the machinery worth it? — baselines</h2>
      <p class="lede">Same folds, same similarity split, same conformal wrapper. The shipped model
      is a random forest on 34 descriptors; the baselines are a logistic regression on those same
      descriptors and a logistic regression on 3-mer counts (order-aware).</p>
      <div class="table-wrap">
      <table class="mini-table wide">
        <thead><tr><th>task</th><th>model</th><th>AUROC</th><th>AUPRC</th></tr></thead>
        <tbody>
        ${["bbb","cpp"].flatMap(o=>{const h=a.tasks[o].similarity,p=[["random forest · 34 descriptors (shipped)",h.auroc_mean,h.auprc_mean]];for(const[d,f]of[["lr_features","logistic regression · 34 descriptors"],["lr_3mer","logistic regression · 3-mer counts"]]){const b=a.tasks[o].baselines?.[d];b&&p.push([f,b.auroc_mean,b.auprc_mean])}return p.map(d=>`<tr><td>${o==="bbb"?"BBB penetration":"Cell penetration"}</td>
              <td>${u(d[0])}</td><td class="num">${d[1].toFixed(3)}</td>
              <td class="num">${d[2].toFixed(3)}</td></tr>`)}).join("")}
        </tbody>
      </table>
      </div>
      <p class="hint">The forest beats a linear model on the same descriptors
      (${a.tasks.bbb.similarity.auroc_mean.toFixed(3)} vs
      ${(a.tasks.bbb.baselines?.lr_features?.auroc_mean??NaN).toFixed(3)} on BBB), so the
      non-linearity is doing real work rather than decorating a linear boundary. The 3-mer model
      is competitive, which is the honest framing: composition and order both carry signal, and
      this feature set only has the former.</p>
      ${a.tasks.bbb.similarity.auroc_vs_cpp_mean!==void 0?`<p class="hint">BBB negatives are a mix of hard (CPP) and easy (random) decoys, and the
      headline AUROC hides it: against CPP negatives the same model scores
      ${a.tasks.bbb.similarity.auroc_vs_cpp_mean.toFixed(3)}, against random negatives
      ${(a.tasks.bbb.similarity.auroc_vs_random_mean??NaN).toFixed(3)} (mixed
      ${a.tasks.bbb.similarity.auroc_mean.toFixed(3)}). Do not read the mixed number as a
      natural-library barrier.</p>`:""}
    </section>

    ${a.order_control?`<section class="card">
      <h2>Order control: the scramble, measured</h2>
      <p class="lede">Take a real BBB shuttle and shuffle its residues
      (<code>${u(a.order_control.parent)}</code> →
      <code>${u(a.order_control.scramble)}</code> — the same 19 residues). Any model that only
      reads composition must give the identical answer; only a model that reads order can tell
      them apart.</p>
      <div class="table-wrap">
      <table class="mini-table wide">
        <thead><tr><th>model</th><th>parent</th><th>scramble</th><th>|ΔP|</th></tr></thead>
        <tbody>${Object.entries(a.order_control.models).map(([o,h])=>`<tr><td>${u(o==="rf"?"shipped random forest":o==="lr_features"?"logistic regression · 34 descriptors":"logistic regression · 3-mer counts")}</td>
          <td class="num">${h.parent.toFixed(3)}</td>
          <td class="num">${h.scramble.toFixed(3)}</td>
          <td class="num">${h.abs_gap.toFixed(3)}</td></tr>`).join("")}</tbody>
      </table>
      </div>
      <p class="hint">${u(a.order_control.note)} This is the sharpest statement of the model's
      limit: a scrambled Angiopep-2 is not a BBB shuttle, and the shipped features cannot express
      that. The fix is an order-aware or learned-embedding feature path, not more trees.</p>
    </section>`:""}

    <section class="card warn">
      <h2>What this benchmark does not prove</h2>
      <ul class="reasons">
        <li>We could not re-run B3Pred itself under similarity splits: its standalone package
        references a server-side model file (<code>/gpsr/models/bp3red/RF_model</code>) that is
        not distributed. Published numbers are quoted, not reproduced.</li>
        <li>The BBB positives are ${a.dataset.bbb_pos} experimentally validated peptides.
        Every metric here is a small-sample estimate; the similarity split makes it honest,
        not precise.</li>
        <li>Negatives are random and CPP peptides, so the model separates BBB-penetrating
        peptides from <em>those</em> decoys — not from a natural, drug-like non-BBB library.</li>
        <li>Cell penetration is trained on CPP-vs-random data as a membrane-permeability
        proxy; it says nothing about tumour selectivity.</li>
        <li>No learned-embedding baseline (e.g. ESM) is included. It needs torch and ~600 MB of
        weights and cannot run in a browser-only tool; the 3-mer baseline answers the same
        question — "is the descriptor machinery worth it?" — at a fraction of that cost.</li>
      </ul>
    </section>`,ne(t.querySelector("#bars-a"),{groups:n,seriesNames:["B3Pred published RF","ours (their split)","ours (novel-only)"],max:1,yLabel:"AUROC",annotations:r}),ne(t.querySelector("#bars-b"),{groups:m,seriesNames:["AUROC","AUPRC"],max:1,yLabel:"score"}),se(t.querySelector("#line-abstain"),{series:[{name:"BBB · similarity split",points:c("bbb","similarity")},{name:"BBB · random split",points:c("bbb","random"),dash:!0},{name:"Cell · similarity split",points:c("cpp","similarity")}],xLabel:"abstention rate",yLabel:"error among answered",xDomain:[0,.7],yDomain:[0,.3]}),se(t.querySelector("#line-reliability"),{series:[{name:"BBB · similarity",points:a.tasks.bbb.similarity.reliability.map(o=>({x:o.predicted,y:o.observed}))},{name:"Cell · similarity",points:a.tasks.cpp.similarity.reliability.map(o=>({x:o.predicted,y:o.observed})),dash:!0}],xLabel:"predicted probability",yLabel:"observed frequency",xDomain:[0,1],yDomain:[0,1],diagonal:!0})}function Je(t,e){const a=B(e.ood_sim);t.innerHTML=`
    <section class="card">
      <h2>What this is</h2>
      <p class="lede">A batch triage tool for glioblastoma peptide candidates. It scores four
      objectives, draws the Pareto trade-offs, and — the part served by no existing BBB
      server — tells you when it does not know, with a calibrated reason.</p>
      <div class="grid-2">
        <div class="card inner">
          <h3>Trained and calibrated</h3>
          <ul class="reasons">
            <li><strong>BBB penetration</strong> — random forest on B3Pdb/B3Pred data,
            conformal prediction sets with a marginal coverage guarantee.</li>
            <li><strong>Cell penetration</strong> — same machinery on CPP-vs-random data;
            a membrane-permeability proxy, not a tumour-selectivity claim.</li>
          </ul>
        </div>
        <div class="card inner">
          <h3>Transparent heuristic indexes</h3>
          <ul class="reasons">
            <li><strong>Immune compatibility</strong> — MHC-II anchor density and class-II
            length window, with every component shown.</li>
            <li><strong>Manufacturability</strong> — charge, aggregation runs and known
            liability motifs (deamidation, oxidation, isomerisation, free cysteines).</li>
          </ul>
        </div>
      </div>
      <p class="hint">Heuristic indexes are deterministic rules, deliberately not sold as
      predictions. For calibrated MHC-II binding, run the sequences through NetMHCIIpan and
      override.</p>
    </section>

    <section class="card">
      <h2>What a score means — and what it does not</h2>
      <p class="lede">Read this before quoting any number from this tool.</p>
      <div class="grid-2">
        <div class="card inner">
          <h3>What it means</h3>
          <ul class="reasons">
            <li><strong>P(BBB) is a calibrated ranking score on a specific comparison</strong> —
            validated BBB-penetrating peptides against CPP and random decoys. It is not a
            probability of crossing the human blood–brain barrier.</li>
            <li><strong>A confident call means the two labels separated</strong> on held-out
            peptide clusters at the chosen coverage target — nothing about in vivo efficacy.</li>
            <li><strong>Abstention is a result, not an error.</strong> It says the two labels are
            not separable at that coverage on data this model has seen.</li>
            <li><strong>Similarity and the out-of-domain flag</strong> tell you whether the
            peptide sits inside the region where the estimate is meaningful at all.</li>
          </ul>
        </div>
        <div class="card inner">
          <h3>What it does not mean</h3>
          <ul class="reasons">
            <li>Not measured permeability: no assay, no P<sub>app</sub>, no
            K<sub>p,uu,brain</sub>.</li>
            <li>Not the endothelium, efflux transporters (P-gp, BCRP), plasma-protein binding,
            proteolytic stability, or dose.</li>
            <li>Cell penetration is not tumour selectivity, and neither is brain delivery.</li>
            <li>The immune and manufacturability indexes are deterministic rules, not
            immunogenicity or developability assays.</li>
            <li>A GBM-target hit is a sequence lookup, not evidence of binding.</li>
            <li>A positive here does not accelerate drug discovery — it prioritises which
            peptides are worth testing.</li>
          </ul>
        </div>
      </div>
    </section>

    <section class="card">
      <h2>Conformal prediction (the abstention)</h2>
      <p class="lede">Split conformal with the LAC nonconformity score. For each candidate
      sequence x and label y, score s(x, y) = 1 − p̂(y | x). On a held-out calibration set of
      n peptides, take q as the ⌈(n+1)(1−α)⌉/n empirical quantile of the true-label scores.
      The prediction set is {y : s(x, y) ≤ q}.</p>
      <ul class="reasons">
        <li>One label survives → a call is made (BBB+ / BBB−).</li>
        <li>Both survive → <strong>abstain</strong>: under the coverage target, this peptide
        cannot be separated.</li>
        <li>Coverage is marginal and in-distribution: across many draws, the true label is in
        the set at least 1−α of the time. It is a guarantee on the set, not a promise about
        any single peptide.</li>
      </ul>
      <p class="hint">The α dial in the Predict tab changes q only — the numbers next to it are
      measured on held-out clusters, not asserted. One honest caveat: this is split conformal
      on a single 158-peptide calibration set, so the coverage guarantee is marginal and
      seed-dependent, and the six-point α grid lands on discrete empirical quantiles — that is
      why the measured dial curve is not monotone step to step.</p>
    </section>

    <section class="card">
      <h2>Why the tool can refuse: out-of-domain</h2>
      <p class="lede">A conformal guarantee is void for peptides that do not resemble the
      training set. For every input we compute sequence similarity to all training peptides,
      sim(a,b) = 2·LCS(a,b) / (|a|+|b|), and compare the maximum against a gate of
      <strong>${a}</strong>. Below that — or when the length falls outside the training
      range of 5–30 aa — the peptide is flagged out-of-domain and the conformal verdict is
      withheld: the interface says so instead of printing a confident number.</p>
      <p class="hint">The gate is not a guessed constant. It is the 10th percentile of the
      similarity from each held-out test peptide to the training pool (0.43 here, floored at
      0.40), so by construction it flags the bottom decile of peptides that already resemble
      the training data. A fixed 0.40 was tried first and measured: it flagged only 0.6% of
      held-out peptides, i.e. it never fired on realistic input. The number and its failure
      rate are both exported in <code>benchmark.json</code>.</p>
    </section>

    <section class="card">
      <h2>The models, exactly</h2>
      <p>Random forest, 200 trees, max depth 6, min 3 samples per leaf, class-balanced, on
      ${z.length} transparent descriptors per peptide:</p>
      <p class="mono wrap-text">${z.join(" · ")}</p>
      <p class="hint">No 9,000-feature Pfeature kitchen sink, no learned embeddings. The point
      of this tool is not raw accuracy — it is knowing when accuracy is not there. Full
      disclosure: our BBB forest on the similarity split reaches AUROC ≈ 0.82, below B3Pred's
      published 0.88–0.90 on random splits, and part of that gap is the split protocol itself
      (Validation tab).</p>
    </section>

    <section class="card">
      <h2>Heuristic definitions (so you can disagree precisely)</h2>
      <div class="grid-2">
        <div class="card inner">
          <h3>Immune compatibility (higher = safer)</h3>
          <p>Start at 100, subtract:</p>
          <ul class="reasons">
            <li><strong>75 ×</strong> mean MHC-II anchor density — over every 9-mer frame,
            the fraction of P1/P4/P6/P9 positions occupied by a hydrophobic residue
            (A/V/I/L/M/F/W/Y).</li>
            <li><strong>25 ×</strong> length factor — 1 inside 13–25 aa (optimal class-II
            window), 0.5 for 9–12 or 26–30, else 0.15.</li>
          </ul>
        </div>
        <div class="card inner">
          <h3>Manufacturability (higher = easier)</h3>
          <p>50 × charge score (favorable |net charge| 2–6 at pH 7.4, falling off outside)
          + 50 × aggregation score (longest hydrophobic run ≤ 4 good, ≥ 8 bad), minus
          penalties: −8 per deamidation motif (N[GSTHA]), −4 per Met, −6 per Trp, −8 per
          free cysteine, −6 per Asp isomerisation motif (D[GST]), −6 for N-terminal Gln.
          Each penalty caps at 4 occurrences.</p>
        </div>
      </div>
    </section>

    <section class="card">
      <h2>GBM-target lookup (a lookup, not a prediction)</h2>
      <p>Glioblastoma is the point of this tool, and a BBB-penetrating peptide is only useful if
      it also goes somewhere. The GBM column is an explicit, cited substring lookup against
      published GBM-relevant binder sequences and adhesion motifs — nothing is inferred from the
      model. A hit means the sequence <em>contains</em> a published binder's sequence; affinity,
      selectivity, stability and actual brain delivery are not implied.</p>
      <div class="table-wrap">
      <table class="mini-table wide">
        <thead><tr><th>entry</th><th>sequence match</th><th>target</th><th>citation</th></tr></thead>
        <tbody>${he.map(i=>`<tr><td>${u(i.name)}</td><td class="mono">${u(i.pattern)}</td>
          <td>${u(i.target)}</td><td class="hint">${u(i.cite)}</td></tr>`).join("")}</tbody>
      </table>
      </div>
      <p class="hint">Deliberately excluded, because they have no short linear motif to match:
      ${pe.map(i=>`<strong>${u(i.target)}</strong> (${u(i.why)})`).join("; ")}.
      A missing hit therefore means "not in this lookup", not "not a GBM target".</p>
    </section>

    <section class="card">
      <h2>Benchmark protocol</h2>
      <ul class="reasons">
        <li><strong>Random split</strong> — stratified 60/20/20, the protocol the field
        reports.</li>
        <li><strong>Similarity split</strong> — every peptide clustered at ≥80% mutual
        LCS-similarity (union-find over all pairs); the classifier never sees a cluster in
        training that appears in test. Calibration uses separate held-out clusters.</li>
        <li><strong>B3Pred audit</strong> — our model trained on B3Pred's published train
        partitions, evaluated on their published validation partitions; then the similarity
        of every validation peptide to the training pool is measured.</li>
        <li>All numbers: 5 seeds, mean ± sd, reproduced by <code>python3 scripts/train.py</code>.</li>
        <li><strong>Baselines</strong> — same folds, same conformal wrapper: logistic
        regression on the same 34 descriptors, and logistic regression on 3-mer counts
        (order-aware, vocabulary fitted on the training fold only). The forests beat both on the
        similarity split, which is the check that the extra machinery is earning its keep.</li>
        <li><strong>Per-negative-type AUROC</strong> — the BBB negatives are a mix of CPPs (hard)
        and random peptides (easy). The headline AUROC averages over both, so it is reported
        broken out as well: a strong number against random decoys is not a natural-library
        barrier.</li>
        <li><strong>Coverage CI</strong> — a Wilson interval on the empirical test coverage, so a
        gap below the nominal target can be judged against sampling noise instead of read as a
        bug.</li>
        <li><strong>Order control</strong> — a real shuttle and an anagram of it are pushed
        through every model; the permutation-invariant ones must agree exactly.</li>
      </ul>
    </section>

    <section class="card">
      <h2>Known failure cases (read before trusting anything)</h2>
      <ul class="reasons">
        <li><strong>Non-canonical residues</strong> (D-amino acids, N-methylation,
        β-amino acids, cyclisation, PEGylation) are invisible to every descriptor here — the
        input validation rejects them rather than pretending.</li>
        <li><strong>Lengths outside 5–30 aa</strong> extrapolate; the UI flags them as
        out-of-domain.</li>
        <li><strong>Scrambled sequences beat composition-only features.</strong> 33 of the 34
        descriptors are permutation-invariant (the 20 residue fractions, length, net charge,
        GRAVY, MW and the crude class fractions); the only order-sensitive one is the longest
        hydrophobic run, and a shuffle usually keeps that too. The example batch includes a
        scrambled Angiopep-2 control and the model gives it the <em>identical</em> probability
        (measured |ΔP| = 0.000 — see the order control on the Validation tab), while a 3-mer
        baseline separates the pair by 0.22. Position-aware features are the fix; this tool does
        not have them, so treat high probabilities on shuffled or low-complexity sequences as
        uninformative.</li>
        <li><strong>Membrane-active CPPs</strong> that are not BBB shuttles can score high on
        cell penetration and moderate on BBB — judged against CPP/random negatives, not
        against the BBB endothelium in vivo.</li>
        <li><strong>Aggregation-prone leads</strong> may pass the manufacturability index:
        the hydrophobic-run rule is a screening heuristic, not a Zyggregator-class predictor.</li>
        <li><strong>Immunogenicity index ≠ immunogenicity.</strong> HLA-DRB1 alleles differ;
        a two-anchor rule cannot capture allele-specific binding.</li>
        <li><strong>Small positives (265 after removing sequences duplicated in the CPP
        negatives; 269 in the source file).</strong> The similarity-split BBB number is the one
        that matters, and it is 0.82 AUROC, not 0.95.</li>
        <li><strong>Short peptides get similarity for free.</strong> 2·LCS/(|a|+|b|) is large
        for short sequences (a 3-mer shares its whole self with any 5-mer), so the
        out-of-domain gate is weakest exactly where the feature set is also weakest.</li>
      </ul>
    </section>

    <section class="card">
      <h2>Data & citations</h2>
      <ul class="reasons">
        <li>Positive/negative data: B3Pred web server datasets (B3Pdb-derived), downloaded
        ${new Date().toISOString().slice(0,10)}. Kumar V., Patiyal S., Dhall A., Sharma N.,
        Raghava G.P.S. <em>B3Pred: A Random-Forest-Based Method for Predicting and Designing
        Blood–Brain Barrier Penetrating Peptides.</em> Pharmaceutics 2021, 13(8):1237.</li>
        <li>Conformal prediction: Lei J. et al. <em>Distribution-Free Predictive Inference
        for Regression.</em> JASA 2018; Angelopoulos & Bates, <em>A Gentle Introduction to
        Conformal Prediction and Distribution-Free Uncertainty Quantification</em></li>
        <li>Kyte–Doolittle hydropathy (1982); LCS-based sequence identity as defined above.
        Ramírez et al. for the minimal-prefix notation this deliberately avoids.</li>
      </ul>
      <p class="hint">Not built yet, not claimed: calibrated immunogenicity (Needs
      NetMHCIIpan), design mode with single-residue edits, full-data refit. Each is a clear
      next step on top of this scaffold, none is faked here.</p>
    </section>`}const Ze=[{id:"predict",label:"Predict",blurb:"batch triage with abstention"},{id:"pareto",label:"Pareto",blurb:"which trade-offs are real"},{id:"validation",label:"Validation",blurb:"does it generalise"},{id:"methods",label:"Methods",blurb:"how it works, and its limits"}],T=document.getElementById("root");function et(){T.innerHTML='<div class="splash">loading model…</div>';const t="./";fetch(`${t}model.json`).then(e=>{if(!e.ok)throw new Error(`model.json: HTTP ${e.status}`);return e.json()}).then(e=>tt(e)).catch(e=>{T.innerHTML=`<div class="splash"><h2>Could not load the model</h2>
        <p class="hint">${u(String(e))}</p>
        <p class="hint">The app needs <code>public/model.json</code>, generated by
        <code>python3 scripts/train.py</code>. If it is missing from the deployed build,
        the model file was not committed or the base path is wrong.</p></div>`}),fetch(`${t}benchmark.json`).then(e=>e.ok?e.json():null).then(e=>{if(e){if(!x){be=e;return}x.benchmark=e,ge(),M==="validation"&&A()}}).catch(()=>{})}let x,M="predict",be=null,O;function tt(t){x={model:t,benchmark:be,input:I,alphaIndex:j.indexOf(.1),results:[],selected:null,picked:null,paretoX:"bbb",paretoY:"cpp",showGhosts:!0},T.innerHTML=`
    <header class="top">
      <div class="title-row">
        <div>
          <h1>Peptide triage for glioblastoma</h1>
          <p class="sub">Blood–brain-barrier and cell penetration, with a conformal
          <em>I don't know</em>. Every number below is either measured on held-out data or a
          named heuristic — nothing is a black box.</p>
        </div>
      </div>
      <nav class="tabs" role="tablist">
        ${Ze.map(a=>`<button role="tab" id="tab-${a.id}" aria-controls="panel" data-tab="${a.id}" class="tab">
              ${a.label}<small>${a.blurb}</small>
            </button>`).join("")}
      </nav>
      <div class="dial-bar">
        <label class="dial">
          <span>coverage target <strong id="dial-coverage">90%</strong>
            <span class="hint">(α = <span id="dial-alpha">0.10</span>)</span></span>
          <input type="range" id="dial" min="0" max="${j.length-1}" step="1"
            value="${x.alphaIndex}" aria-label="conformal coverage target">
          <span class="dial-scale"><i>α = 0.40 · looser</i><i>α = 0.05 · tighter</i></span>
        </label>
        <p class="hint" id="dial-note">Loading measured dial behaviour…</p>
      </div>
    </header>
    <main id="panel" role="tabpanel"></main>
    <footer class="bottom">
      Trained on B3Pdb/B3Pred peptides (Kumar V. et al., Pharmaceutics 2021, 13(8):1237).
      Heuristic indexes are rules, not models. Not a clinical tool — a screening aid that
      shows its reasoning and refuses when it should.
    </footer>`;const e={setInput(a){x.input=a},analyze(){x.results=N(x.model,X(x.input),x.alphaIndex),x.selected=null,x.picked=null,A()},setAlpha(a){x.alphaIndex=a,x.results.length&&(x.results=N(x.model,X(x.input),a)),A()},select(a){x.selected=a,A()},pick(a){x.picked=a,A()},setAxis(a,i){a==="x"?x.paretoX=i:x.paretoY=i,A()},toggleGhosts(){x.showGhosts=!x.showGhosts,A()},reRender:A,goTab(a){M=a,A()}};T.querySelectorAll("[data-tab]").forEach(a=>a.addEventListener("click",()=>{M=a.dataset.tab,A()})),O=e,T.querySelector("#dial").addEventListener("input",a=>e.setAlpha(Number(a.target.value))),document.addEventListener("keydown",a=>{(a.metaKey||a.ctrlKey)&&a.key==="Enter"&&(a.preventDefault(),e.analyze())}),e.analyze()}function ge(){const t=j[x.alphaIndex],e=document.getElementById("dial-coverage"),a=document.getElementById("dial-alpha"),i=document.getElementById("dial-note");if(e&&(e.textContent=B(1-t)),a&&(a.textContent=t.toFixed(2)),!i)return;const s=x.benchmark;if(!s){i.textContent="Benchmark still loading — the dial still works; the measured numbers arrive with benchmark.json.";return}const n=d=>d.reduce((f,b)=>Math.abs(b.alpha-t)<Math.abs(f.alpha-t)?b:f),r=s.tasks.bbb.similarity.alpha_sweep,m=s.tasks.cpp.similarity.alpha_sweep,c=n(r),o=n(m),h=r[0],p=r[r.length-1];i.textContent=`measured on held-out clusters at α=${c.alpha}: BBB leaves ${B(c.abstain_rate)} unanswered and is wrong on ${B(c.error_rate)} of the answers it gives; cell penetration leaves ${B(o.abstain_rate)} unanswered and is wrong on ${B(o.error_rate)}. Across the whole dial (α ${h.alpha} → ${p.alpha}) BBB goes from ${B(h.abstain_rate)} to ${B(p.abstain_rate)} refused and from ${B(h.error_rate)} to ${B(p.error_rate)} wrong — more refusals, fewer mistakes. Individual steps are not monotone: each α lands on a different empirical quantile of a small calibration set.`}function A(){const t=document.getElementById("panel");t&&(ge(),document.querySelectorAll("[data-tab]").forEach(e=>{const a=e.dataset.tab===M;e.classList.toggle("active",a),e.setAttribute("aria-selected",String(a))}),M==="predict"?Oe(t,x,O):M==="pareto"?ze(t,x,O):M==="validation"?Xe(t,x):Je(t,x.model))}et();
