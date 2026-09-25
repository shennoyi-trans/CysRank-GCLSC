import unittest
import torch
from src import core, ranking, pipeline


class RankingTests(unittest.TestCase):
    def test_duplicate_conformers_do_not_increase_weight(self):
        head=torch.nn.Linear(2,1)
        x=torch.tensor([[1.,0.],[0.,1.]])
        y=torch.tensor([1.,0.]);w=torch.ones(2)
        pair={'preferred':x[:1], 'other':x[1:]}
        many={k:v.repeat(100,1) for k,v in pair.items()}
        a=core.joint_loss(head,x,y,w,[0,1],pair)[0]
        b=core.joint_loss(head,x,y,w,[0,1],many)[0]
        self.assertAlmostEqual(a.item(),b.item(),places=6)

    def test_pair_loss_direction(self):
        head=torch.nn.Linear(2,1)
        with torch.no_grad():head.weight.copy_(torch.tensor([[1.,-1.]]));head.bias.zero_()
        x=torch.eye(2);y=torch.tensor([1.,0.]);w=torch.ones(2)
        good=core.joint_loss(head,x,y,w,[0,1],{'preferred':x[:1],'other':x[1:]})[2]
        bad=core.joint_loss(head,x,y,w,[0,1],{'preferred':x[1:],'other':x[:1]})[2]
        self.assertLess(good,bad)

    def test_sst_has_no_binary_labels_and_shared_group(self):
        rows=pipeline.prepare(core.ROOT/'data/raw',core.ROOT/'data/provenance/parent_sequences.jsonl')
        event=ranking.load(core.ROOT/'data/ranking/sst_ordering.json',rows)
        self.assertEqual(event['parent_group'],'old3_positive_05')
        self.assertEqual(len(event['samples']),100)
        self.assertEqual(sum(sum(r['label_mask']) for r in rows),22)
        self.assertTrue(all('ca_label' not in s['sample'] for s in event['samples']))
