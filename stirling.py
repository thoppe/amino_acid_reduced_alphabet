from copy import deepcopy
import itertools

def stirling_set(base, items=list(), index=0, 
                 early_break=None,
                 early_kwargs={}):

    if early_break != None:
        if early_break(items,
                       base=base,index=index,
                       **early_kwargs): 
            raise StopIteration

    # Base case
    if index==len(base):
        yield items
        raise StopIteration

    v = base[index]       

    for k,block in enumerate(items):
        new_items = deepcopy(items)
        new_items[k].append(v)
        for sol in stirling_set(base, new_items, index+1, 
                                early_break, early_kwargs):
            yield sol

    # Add a new set too
    new_items = deepcopy(items)
    new_items.append( [v,] )
    for sol in stirling_set(base, new_items, index+1, 
                            early_break,early_kwargs):
        yield sol


if __name__ == "__main__":

    def max_set_size(A,**kwargs):
        k = 2
        return len(A) > k
    
       
    A = range(4)
    for x in stirling_set(A, [], 
                          early_break=max_set_size,
                          early_kwargs={}):
        print x
