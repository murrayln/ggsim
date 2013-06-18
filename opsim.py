
import random
from UserList import UserList
from collections import OrderedDict
import copy
import scipy 
"""
'group': index of gene group
'pos': index of gene group
For example:
[[4,2,6], [3, 1, 4,7], [8,2]]
    
group 1 [3,1,4,7]
pos 3 in group 1 = 7

"""

class GeneGroups(UserList):
    def __init__(self,inlist=[],lastmut='null'):
        self.data = inlist
        self.lastmut = lastmut
    def duplicate(self,group, pos):
        """ Duplicate a gene in a group"""
        newgroup = self.data[group][:pos] + [self.data[group][pos]] + \
                   self.data[group][pos:]
        retval = self.data[:group] + [newgroup] + self.data[group+1:]
        return GeneGroups(retval,'duplicate')
    def delete(self, group, pos):
        """Delete a gene in a group. If group empty, delete group"""
        retval = copy.deepcopy(self.data)
        del retval[group][pos]
        if len(retval[group]) == 0:
            del retval[group]
        return GeneGroups(retval,'delete')
    def revstrand(self, group, pos):
        """ Move gene to reverse strand, within group"""
        retval = copy.deepcopy(self.data)
        retval[group][pos] = -retval[group][pos]
        # TODO: Need migrate the gene from the group
        # Doesn't make sense to leave it in the group
        return GeneGroups(retval,'revstrand')
    def migrate(self, group, pos):
        """Split a group in two"""
        tmpval = copy.deepcopy(self.data)
        newgroups = [tmpval[group][:pos], tmpval[group][pos:]]
        retval = tmpval[:group] + newgroups + tmpval[group+1:]
        return GeneGroups(retval,'migrate')
    def join(self, group1, group2):
        """Join two gene groups together. Opposite of migrate"""
        # verify group1 is the lower group index
        if group1 > group2:
            group1, group2 = group2, group1
        retval = copy.deepcopy(self.data)
        retval[group1] += retval[group2]
        del retval[group2]
        return GeneGroups(retval,'join')
    def swap(self, group,pos):
        """ Swap two genes in a group """
        retval = copy.deepcopy(self.data)
        if len(retval[group]) > 1:
            exlist = range(len(self.data[group]))
            del exlist[pos]
            pos2 = random.choice(exlist)
            retval[group][pos], retval[group][pos2] = \
            retval[group][pos2], retval[group][pos]
        return GeneGroups(retval,'swap')
    def __repr__(self):
        return str(self.data)
    def __eq__(self, other):
        if other.data == self.data:
            return True
        else:
            return False
    def __ne__(self, other):
        if other.data != self.data:
            return True
        else:
            return False
    def __str__(self):
        return str(self.data)
    def __len__(self):
        return len(self.data)


probmut_lut = OrderedDict([
                ((0,0.005), 'duplicate'),
                ((0.005,0.010), 'delete'),
                ((0.010,0.020), 'revstrand'),
                ((0.020,0.030), 'migrate'),
                ((0.030,0.080), 'swap'),
                ((0.080,0.090), 'join')
                ])

# Fitness function.... not sure. Price Equation?
# For now, I'll try something more like the Weasel program
fitness_lut = {'join': 3,
               'duplicate': 2,
               'delete': -1,
               'revstrand': -2,
               'migrate': -1,
               'swap': 0,
               'null': 0}

def score_brood(brood):
    """ Weighted choice 

     Give a fitness score to each brood member
     Fitness function assumes the following:
     1. joins increase fitness
     2. Gene duplication increases fitness
     3. Deletion decreases fitness
     4. revstrand decreases fitness
     5. Migration decrease fitness
     6. Swaps are neutral
     7. Nulls are neutral
     """

    sorted_brood = []
    for i in brood:
        try:
            sorted_brood.append( (fitness_lut[brood[i].lastmut], 
                              brood[i]) )
        except AttributeError:
            print "AttributeError"
            print i,brood[i]
            raise
    sorted_brood.sort()
    sorted_brood.reverse()
    scored_brood = OrderedDict([ (i,x[1]) for i,x in
                               enumerate(sorted_brood) ])
    return scored_brood
    
def write_scored_brood(generation, scored_brood):
    outhandle = open("brood.%d.csv" % generation, "w") 
    for i in scored_brood:
        lastmut = scored_brood[i].lastmut
        gg = scored_brood[i]


        outhandle.write("%d\t%s\t%s\n" % (i,
                        lastmut,
                        str(gg))
                        )
    outhandle.close()

def read_scored_brood(generation):
    scored_brood = []
    inhandle = open("brood.%d.csv" % generation) 
    for inline in inhandle:
        fitness_score, lastmut, str_offsrping = inline.strip().split('\t')
        offspring = [int(x) for x in str_offspring[1:-1].split(',')]
        scored_brood.append( (fitness_score, lastmut, offspring) )
    return scored_brood
        

def gg_son(gg):
    rnum = random.random()
    event = None
    gg_son = copy.deepcopy(gg)
    for low, high in probmut_lut:
        if rnum >= low and rnum < high:
            event = probmut_lut[(low,high)]
    if event:
        if event == 'join':
            if len(gg) == 1:
                return gg_son
            elif len(gg) == 2:
                group1 = 0; group2 = 1
            elif len(gg) > 2:
                group1 = 0; group2 = 0
                while group1 == group2:
                    group1 = random.choice(range(len(gg)))
                    group2 = random.choice(range(len(gg)))
            print "join", group1, group2
            gg_son = gg.join(group1, group2)
        elif event == 'migrate':
            group = random.choice(range(len(gg)))
            if len(gg[group]) > 1:
                print event, "group:", group, 
                pos = random.choice( range(1,len(gg[group])) )
                print "pos",pos
                gg_son = getattr(gg,event)(group,pos)
        else:
            group = random.choice(range(len(gg)))
            pos = random.choice(range(len(gg[group])))
            print event, "group:", group, "pos:", pos
            gg_son = getattr(gg,event)(group,pos)
    return gg_son
        
def gg_brood(gg,broodsize=1000):
    brood = OrderedDict([])
    for i in range(broodsize):
#        brood.append(gg_son(gg))
        brood[i] = gg_son(gg)
        #print i,brood[i]
    return brood

def multi_gg_brood(in_brood, generation=0, broodsize=1000,breedfrac=0.05):
    next_brood = OrderedDict([])
    culled_broods = []
    breedsize = int(breedfrac * broodsize)
    for gg in in_brood.values():
#        print type(gg); break
        brood = gg_brood(gg, broodsize)    
        # take only the breeding fractions from each brood
        culled_broods.append(score_brood(brood).values()[:breedsize])
    # treat all breeding ggs in this generation as a brood
    i = 0
    
    for brood in culled_broods:
        for gg in brood:
            next_brood[i] = gg
            i += 1
    next_brood = score_brood(next_brood)
    return OrderedDict(next_brood.items()[:broodsize])
    

def run_generations(gg,broodsize=1000,ngen=100,breedfrac=0.05):
    breedsize = int(breedfrac * broodsize)
    brood = gg_brood(gg,broodsize)
    next_brood = score_brood(brood) 
    write_scored_brood(1,next_brood)
    for generation in range(2,ngen+1):
        print "*********************"
        print "generation", generation
        print "*********************"
        next_brood = multi_gg_brood(next_brood, generation)
        write_scored_brood(generation,next_brood)
        
        
