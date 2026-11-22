"""Data and exact Euclidean scoring for hierarchy-aware constellations."""
from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path

import networkx as nx
import numpy as np
import pandas as pd
import torch

ROOT=Path(__file__).resolve().parents[1]
SOURCE=ROOT/'results/hierarchy_global_reanalysis/final_132'
VERSION=int(os.environ.get('CONSTELLATION_VERSION','1'))
assert VERSION in (1,2,3)
OUT=SOURCE/f'constellation_discovery_v{VERSION}'
CONFIG=dict(version=VERSION,inner_seed=20261029,validation_studies=87,graph_seed=41,
            graph_resolution=1.0,top_partners=5,min_discovery_donors=20,
            min_discovery_studies=5,min_validation_donors=10,min_validation_studies=3,
            min_members=2,refinement_improvement=1e-4,
            metric='Arithmetic mean pairwise Euclidean distance on 2,000 training-selected genes',
            pair_score='1 - partner-shifted target distance / unshifted target distance on identical maps',
            graph='Union of discovery top-five edges with positive gain in both directions; mean directional gain; Louvain resolution 1',
            optimization='Minimize equal-member mean leave-one-star-out distance/raw-distance ratio on fixed internal-validation maps',
            hierarchy='Ancestors, descendants and other overlapping types allowed; exact target star excluded from its own prediction',
            coverage='Complete constellation profiles on a fixed seed-specific training benchmark; additions cannot change benchmark maps',
            repeated_views='All maps scored primarily; donor counts gate support and one-view-per-donor and within-study sensitivity are reported',
            literature={'Louvain':'https://arxiv.org/abs/0803.0476',
                        'implementation':'https://networkx.org/documentation/stable/reference/algorithms/generated/networkx.algorithms.community.louvain.louvain_communities.html'})
if VERSION>=2:
    CONFIG['hierarchy']='Exclude every ancestor/descendant pair using transitive full-ontology paths; other overlaps remain eligible; exact target excluded'
    CONFIG['test_status']='Exploratory rerun after examining version-1 results on the same 88 test studies'
if VERSION==3:
    CONFIG['hierarchy']='Exclude ancestors/descendants and pairs with pooled training actual-cell IoU > 0.33; exact target excluded'
    CONFIG['max_cell_iou']=0.33
    CONFIG['iou_scope']='Ratio of summed actual intersection/union counts across common complete outer-training maps; all 348 training studies; no test counts used for eligibility'
    CONFIG['test_status']='Exploratory rerun requested after inspecting versions 1 and 2 on the same 88 test studies'


def digest(path):
    with Path(path).open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()


def write_json(path,value):
    Path(path).write_text(json.dumps(value,indent=2,allow_nan=False)+'\n',encoding='utf-8')


def distance_matrix(values):
    """Double-precision distances; direct subtraction repairs cancellation."""
    x=values.to(dtype=torch.float64)
    if len(x)>0 and torch.equal(x,x[:1].expand_as(x)):
        return torch.zeros((len(x),len(x)),dtype=torch.float64,device=x.device)
    x=x-x.mean(dim=0,keepdim=True)
    norm=(x*x).sum(dim=1)
    d=x@x.T
    d.mul_(-2).add_(norm[:,None]).add_(norm[None,:]).clamp_(min=0)
    d.fill_diagonal_(0)
    positions=torch.nonzero(d<1e-8,as_tuple=False)
    positions=positions[positions[:,0]!=positions[:,1]]
    for start in range(0,len(positions),4096):
        p=positions[start:start+4096]
        d[p[:,0],p[:,1]]=((x[p[:,0]]-x[p[:,1]])**2).sum(dim=1)
    d.sqrt_()
    return d


def mean_distances(values):
    """Input is [type,map,gene] or [map,gene]; no pair subsampling."""
    if values.ndim==2:
        n=len(values)
        if n<2:return float('nan')
        return float(distance_matrix(values).sum().item()/(n*(n-1)))
    return np.asarray([mean_distances(v) for v in values],dtype=np.float64)


def ratios(distance,raw):
    return np.divide(distance,raw,out=np.full_like(np.asarray(raw,dtype=float),np.nan),where=np.asarray(raw)>1e-12)


class Data:
    def __init__(self):
        self.meta=pd.read_csv(SOURCE/'star_rows.csv',low_memory=False)
        self.maps=pd.read_csv(SOURCE/'maps.csv',low_memory=False)
        assert len(self.meta)==875669 and len(self.maps)==16524
        assert np.array_equal(self.meta.star_row,np.arange(len(self.meta)))
        self.labels=sorted(self.meta.parent_id.unique());assert len(self.labels)==132
        self.codes=pd.Categorical(self.meta.parent_id,categories=self.labels).codes
        types=pd.read_csv(SOURCE/'reference_map_types.csv').sort_values('reference_row')
        assert types.parent_id.tolist()==self.labels
        self.names=types['name'].tolist()
        self.reference=np.load(SOURCE/'training_type_reference.npy').astype(np.float64)
        assert self.reference.shape==(132,2000) and np.isfinite(self.reference).all()
        self.raw=np.memmap(SOURCE/'raw_stars.float32',mode='r',dtype='float32',shape=(len(self.meta),2000))
        complete=self.meta.complete_gene_axis.to_numpy(dtype=bool)
        self.lookup=np.full((len(self.maps),132),-1,dtype=np.int64)
        rows=np.flatnonzero(complete)
        self.lookup[self.meta.iloc[rows].map_index,self.codes[rows]]=rows
        self.study_code=pd.Categorical(self.maps.study_id).codes
        self.donor_code=pd.Categorical(self.maps.donor_uid).codes
        frozen=pd.read_csv(ROOT/'results/pooled_star_translation/maps.csv',usecols=['study_id','split']).drop_duplicates()
        assert frozen.study_id.is_unique
        assert self.maps.split.eq(self.maps.study_id.map(frozen.set_index('study_id').split)).all()
        assert self.maps.groupby('donor_uid').split.nunique().le(1).all()
        train_studies=sorted(self.maps.loc[self.maps.split.eq('train'),'study_id'].unique())
        assert len(train_studies)==348
        np.random.default_rng(CONFIG['inner_seed']).shuffle(train_studies)
        validation=set(train_studies[:CONFIG['validation_studies']])
        self.maps['inner_split']=np.where(self.maps.split.ne('train'),'outer_test',
                                         np.where(self.maps.study_id.isin(validation),'validation','discovery'))
        self.fold_maps={fold:np.flatnonzero(self.maps.inner_split.eq(fold)) for fold in ['discovery','validation','outer_test']}
        assert self.maps.iloc[self.fold_maps['outer_test']].study_id.nunique()==88
        self.train_rows=np.flatnonzero(complete & self.meta.split.eq('train').to_numpy())
        self.global_to_train=np.full(len(self.meta),-1,dtype=np.int64)
        self.global_to_train[self.train_rows]=np.arange(len(self.train_rows))
        self.train_lookup=np.full_like(self.lookup,-1)
        mask=self.lookup>=0
        self.train_lookup[mask]=self.global_to_train[self.lookup[mask]]
        self.global_shift=np.load(SOURCE/'globally_shifted_geometry/map_mean_shifts.npy',mmap_mode='r')
        self.map_k=(self.lookup>=0).sum(axis=1)
        members=pd.read_csv(ROOT/'references/vocabulary/aggregation_membership.tsv',sep='\t')
        membership=members.groupby('aggregate_id').original_label.agg(set).to_dict()
        self.members=[membership[x] for x in self.labels]
        edges=pd.read_csv(ROOT/'references/vocabulary/hierarchy_edges.tsv',sep='\t')
        graph=nx.DiGraph();graph.add_edges_from(zip(edges.parent_id,edges.child_id))
        self.ancestors=[nx.ancestors(graph,x) for x in self.labels]
        self.iou_allowed=None;self.training_cell_iou=None
        if VERSION==3 and (OUT/'cell_iou_admissibility.npz').exists():
            saved=np.load(OUT/'cell_iou_admissibility.npz')
            assert saved['labels'].tolist()==self.labels
            self.iou_allowed=saved['allowed'];self.training_cell_iou=saved['pooled_iou']

    def relation(self,target,anchor):
        a,b=self.members[target],self.members[anchor]
        if a==b:kind='identical_direct_membership'
        elif self.labels[anchor] in self.ancestors[target]:kind='anchor_is_ancestor'
        elif self.labels[target] in self.ancestors[anchor]:kind='anchor_is_descendant'
        elif a&b:kind='overlapping_nonancestor'
        else:kind='disjoint_membership'
        result=dict(relationship=kind,shared_direct_labels=len(a&b),direct_label_jaccard=len(a&b)/len(a|b))
        if self.training_cell_iou is not None:
            result['training_cell_iou']=float(self.training_cell_iou[target,anchor])
        return result

    def ontology_allowed(self,a,b):
        return a!=b and (self.labels[a] not in self.ancestors[b] and self.labels[b] not in self.ancestors[a])

    def allowed(self,a,b):
        if VERSION==1:return a!=b
        if not self.ontology_allowed(a,b):return False
        if VERSION==3:
            assert self.iou_allowed is not None,'Prepare training-cell IoU first'
            return bool(self.iou_allowed[a,b])
        return True

    def valid_group(self,members):
        import itertools
        return all(self.allowed(a,b) for a,b in itertools.combinations(members,2))

    def common(self,members,fold):
        maps=self.fold_maps[fold]
        return maps[(self.lookup[np.ix_(maps,list(members))]>=0).all(axis=1)]

    def support(self,maps):
        return dict(maps=len(maps),donors=len(np.unique(self.donor_code[maps])),studies=len(np.unique(self.study_code[maps])))

    def supported(self,maps,fold):
        count=self.support(maps)
        prefix='discovery' if fold=='discovery' else 'validation'
        return count['donors']>=CONFIG['min_'+prefix+'_donors'] and count['studies']>=CONFIG['min_'+prefix+'_studies']

    def values(self,members,maps,device='cuda:0'):
        rows=self.lookup[np.ix_(maps,list(members))]
        assert (rows>=0).all()
        values=np.asarray(self.raw[rows.ravel()],dtype=np.float64).reshape(len(maps),len(members),2000)
        return torch.as_tensor(values.transpose(1,0,2).copy(),device=device)

    def donor_maps(self,members,maps):
        rows=self.lookup[np.ix_(maps,list(members))]
        support=self.meta.n_cells.to_numpy()[rows].min(axis=1)
        table=self.maps.iloc[maps][['study_id','donor_uid','sample_id']].copy()
        table['support']=support;table['original_map_index']=maps
        table=table.sort_values(['support','sample_id'],ascending=[False,True])
        return np.sort(table.drop_duplicates(['study_id','donor_uid']).original_map_index.to_numpy())


def local_metrics(values):
    """Target-excluded mean shift, with hierarchy overlaps explicitly allowed."""
    k=values.shape[0];assert k>=2
    raw=mean_distances(values)
    residual=(k*values-values.sum(dim=0,keepdim=True))/(k-1)
    local=mean_distances(residual)
    ratio=ratios(local,raw)
    return dict(raw=raw,local=local,ratio=ratio,objective=float(np.nanmean(ratio)))


def group_metrics(data,members,maps,device='cuda:0'):
    values=data.values(members,maps,device)
    local=local_metrics(values)
    raw=local['raw'];k=len(members)
    g=torch.as_tensor(np.asarray(data.global_shift[maps]).copy(),device=device,dtype=torch.float64)
    n=torch.as_tensor(data.map_k[maps],device=device,dtype=torch.float64)
    ref=torch.as_tensor(data.reference[list(members)],device=device)
    global_inclusive=mean_distances(values-g[None,:,:])
    global_leaveout=values-(n[None,:,None]*g[None,:,:]-(values-ref[:,None,:]))/(n[None,:,None]-1)
    global_loo=mean_distances(global_leaveout)
    result=[]
    for i,t in enumerate(members):
        result.append(dict(type_index=int(t),parent_id=data.labels[t],name=data.names[t],
                           raw_distance=float(raw[i]),local_loo_distance=float(local['local'][i]),
                           local_inclusive_distance=float(local['local'][i]*(k-1)/k),
                           global_loo_distance=float(global_loo[i]),global_inclusive_distance=float(global_inclusive[i]),
                           local_reduction=float(1-local['ratio'][i]),
                           local_gain_over_global_loo=float(1-local['local'][i]/global_loo[i]) if global_loo[i]>1e-12 else None))
    del values
    return result
