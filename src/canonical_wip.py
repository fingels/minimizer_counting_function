from src.utils import *
from copy import copy, deepcopy

class EnumerationInstance(object):

    __slots__ = 'minimizer','indeterminate_string', 'kmer', 'reverse_complement','m', 'current_position', 'min_starting_position','k', 'alphabet','remaining_letters'

    def __init__(self,string,starting_position,k,alphabet = {'A','C','G','T'}):
        self.m: int = len(string)
        self.k:int = k
        self.alphabet: set[str] = copy(alphabet)
        self.minimizer: str = copy(string)
        self.min_starting_position: int = starting_position
        self.remaining_letters: int = self.k - self.m

        self.indeterminate_string : list[set[str]] = [copy(alphabet)]*k
        self.kmer: list[str] = ['_']*k
        self.reverse_complement: list[str] = ['_']*k

        for j in range(self.m):
            self.indeterminate_string[starting_position + j] = set(string[j])
            self.kmer[starting_position + j] = string[j]
            self.reverse_complement[k - (starting_position + j) - 1] = dna_reverse(string[j])

        self.current_position = 0

    def __str__(self):
        string = 'indeterminate string: '
        for elt in self.indeterminate_string:
            string+= ''.join(sorted(elt))+'|'
        string = string[:-1]+'\n'
        string += '\t\t\t   k-mer: '+' '.join(self.kmer)+'\n'
        string+= '  reverse complement: '+' '.join(self.reverse_complement)+'\n'
        string+='\tcurrent position: '+' '*self.current_position*2+'^'
        return string

    def __copy__(self):

        new_obj = EnumerationInstance(self.minimizer,self.min_starting_position,self.k,self.alphabet)

        new_obj.indeterminate_string = deepcopy(self.indeterminate_string)
        new_obj.kmer = deepcopy(self.kmer)
        new_obj.reverse_complement = deepcopy(self.reverse_complement)

        new_obj.current_position = self.current_position
        new_obj.remaining_letters = self.remaining_letters

        return new_obj

    def process_state(self):

        a = self.kmer[self.current_position]
        b = self.reverse_complement[self.current_position]

        if a is not '_':
            if b is not '_':
                if a<b:
                    # kmer < reverse_complement : terminal state, we can proceed
                    if self.remaining_letters==0:
                        return [('check m-mers',self)]
                    else:
                        return [('complete k-mer', self)]
                elif a==b:
                    # we continue to parse the kmer
                    self.current_position+=1
                    if self.current_position == self.k:
                        # not canonical since kmer = reverse_complement
                        return [('not canonical', 0)]
                    else:
                        return [('keep', self)]
                else:
                    # not canonical since kmer > reverse_complement
                    return [('not canonical',0)]
            else:
                new_instances = []

                #Firt choice : b = a (i.e. dna_reverse(b)=dna_reverse(a))

                new_obj = copy(self)
                new_obj.reverse_complement[new_obj.current_position] = a
                new_obj.kmer[new_obj.k-new_obj.current_position-1] = dna_reverse(a)
                new_obj.indeterminate_string[new_obj.k-new_obj.current_position-1] = set(dna_reverse(a))
                new_obj.current_position+=1
                new_obj.remaining_letters-=1

                if new_obj.current_position==new_obj.k:
                    new_instances.append(('not canonical',0))
                else:
                    new_instances.append(('keep',new_obj))

                # Other choices: b > a (i.e. dna_reverse(b) < dna_reverse(a))
                candidates = {s for s in self.alphabet if s > a}
                for s in candidates:
                    new_obj = copy(self)
                    new_obj.reverse_complement[new_obj.current_position] = s
                    new_obj.kmer[new_obj.k - new_obj.current_position - 1] = dna_reverse(s)
                    new_obj.indeterminate_string[new_obj.k - new_obj.current_position - 1] = set(dna_reverse(s))
                    new_obj.remaining_letters -= 1

                    if new_obj.remaining_letters == 0:
                        new_instances.append(('check m-mers', new_obj))
                    else:
                        new_instances.append(('complete k-mer', new_obj))

                return new_instances
        else:
            #TODO
            return [('did nothing', self)]

def canonical_enumeration(minimizer,k):

    stack: list[EnumerationInstance] = []
    terminal_states = []

    print('Initialization')
    print('**************')

    for i in range(0, k - len(minimizer) + 1):
        obj = EnumerationInstance(minimizer, i, k)

        stack.append(obj)

        # print(obj)
        # print('----')

    print('Processing stack')
    print('****************')

    states_processed = 0

    while len(stack)>0:

        obj = stack.pop()
        states_processed+=1

        # print(obj)

        list = obj.process_state()

        for flag,obj in list:

            # print(flag)

            if flag=='keep':
                stack.append(obj)
            elif flag=='not canonical':
                pass
            else:
                terminal_states.append((flag,obj))

        # print('------')

    return terminal_states, states_processed

minimizer = 'ACA'
k=5
terminal_states, states_processed = canonical_enumeration(minimizer,k)

print('Processed '+str(states_processed)+' internal states.')
print('Number of terminal states: '+str(len(terminal_states)))

for flag,obj in terminal_states:
    print('\n*** '+flag+' ***')
    print(obj)