# (C) 2018 University of Bristol, Bar-Ilan University. See License.txt

import itertools, time
from collections import defaultdict, deque
from Compiler.exceptions import *
from Compiler.config import *
from Compiler.instructions import *
from Compiler.instructions_base import *
from Compiler.util import *
import Compiler.graph
import Compiler.program
import heapq, itertools
import operator
import sys


class StraightlineAllocator:
    """Allocate variables in a straightline program using n registers.
    It is based on the precondition that every register is only defined once."""
    def __init__(self, n):
        self.alloc = {}
        self.usage = Compiler.program.RegType.create_dict(lambda: 0)
        self.defined = {}
        self.dealloc = set()
        self.n = n

    def alloc_reg(self, reg, free):
        base = reg.vectorbase
        if base in self.alloc:
            # already allocated
            return

        reg_type = reg.reg_type
        size = base.size
        if free[reg_type, size]:
            res = free[reg_type, size].pop()
        else:
            if self.usage[reg_type] < self.n:
                res = self.usage[reg_type]
                self.usage[reg_type] += size
            else:
                raise RegisterOverflowError()
        self.alloc[base] = res

        if base.vector:
            for i,r in enumerate(base.vector):
                r.i = self.alloc[base] + i
        else:
            base.i = self.alloc[base]

    def dealloc_reg(self, reg, inst, free):
        self.dealloc.add(reg)
        base = reg.vectorbase

        if base.vector and not inst.is_vec():
            for i in base.vector:
                if i not in self.dealloc:
                    # not all vector elements ready for deallocation
                    return
        free[reg.reg_type, base.size].add(self.alloc[base])
        if inst.is_vec() and base.vector:
            for i in base.vector:
                self.defined[i] = inst
        else:
            self.defined[reg] = inst

    def process(self, program, alloc_pool):
        for k,i in enumerate(reversed(program)):
            unused_regs = []
            for j in i.get_def():
                if j.vectorbase in self.alloc:
                    if j in self.defined:
                        raise CompilerError("Double write on register %s " \
                                            "assigned by '%s' in %s" % \
                                                (j,i,format_trace(i.caller)))
                else:
                    # unused register
                    self.alloc_reg(j, alloc_pool)
                    unused_regs.append(j)
            if unused_regs and len(unused_regs) == len(i.get_def()):
                # only report if all assigned registers are unused
                print "Register(s) %s never used, assigned by '%s' in %s" % \
                    (unused_regs,i,format_trace(i.caller))

            for j in i.get_used():
                self.alloc_reg(j, alloc_pool)
            for j in i.get_def():
                self.dealloc_reg(j, i, alloc_pool)

            if k % 1000000 == 0 and k > 0:
                print "Allocated registers for %d instructions at" % k, time.asctime()

        # print "Successfully allocated registers"
        # print "modp usage: %d clear, %d secret" % \
        #     (self.usage[Compiler.program.RegType.ClearModp], self.usage[Compiler.program.RegType.SecretModp])
        # print "GF2N usage: %d clear, %d secret" % \
        #     (self.usage[Compiler.program.RegType.ClearGF2N], self.usage[Compiler.program.RegType.SecretGF2N])
        return self.usage


def determine_scope(block, options):
    last_def = defaultdict(lambda: -1)
    used_from_scope = set()

    def find_in_scope(reg, scope):
        if scope is None:
            return False
        elif reg in scope.defined_registers:
            return True
        else:
            return find_in_scope(reg, scope.scope)

    def read(reg, n):
        if last_def[reg] == -1:
            if find_in_scope(reg, block.scope):
                used_from_scope.add(reg)
                reg.can_eliminate = False
            else:
                print 'Warning: read before write at register', reg
                print '\tline %d: %s' % (n, instr)
                print '\tinstruction trace: %s' % format_trace(instr.caller, '\t\t')
                print '\tregister trace: %s' % format_trace(reg.caller, '\t\t')
                if options.stop:
                    sys.exit(1)

    def write(reg, n):
        if last_def[reg] != -1:
            print 'Warning: double write at register', reg
            print '\tline %d: %s' % (n, instr)
            print '\ttrace: %s' % format_trace(instr.caller, '\t\t')
            if options.stop:
                sys.exit(1)
        last_def[reg] = n

    for n,instr in enumerate(block.instructions):
        outputs,inputs = instr.get_def(), instr.get_used()
        for reg in inputs:
            if reg.vector and instr.is_vec():
                for i in reg.vector:
                    read(i, n)
            else:
                read(reg, n)
        for reg in outputs:
            if reg.vector and instr.is_vec():
                for i in reg.vector:
                    write(i, n)
            else:
                write(reg, n)

    block.used_from_scope = used_from_scope
    block.defined_registers = set(last_def.iterkeys())

class Merger:
    def __init__(self, block, options, merge_classes):
        self.block = block
        self.instructions = block.instructions
        self.options = options
        if options.max_parallel_open:
            self.max_parallel_open = int(options.max_parallel_open)
        else:
            self.max_parallel_open = float('inf')
        self.dependency_graph(merge_classes)

    def do_merge(self, merges_iter):
        """ Merge an iterable of nodes in G, returning the number of merged
        instructions and the index of the merged instruction. """
        instructions = self.instructions
        mergecount = 0
        try:
            n = next(merges_iter)
        except StopIteration:
            return mergecount, None

        def expand_vector_args(inst):
            new_args = []
            for arg in inst.args:
                if inst.is_vec():
                    arg.create_vector_elements()
                    for reg in arg:
                        new_args.append(reg)
                else:
                    new_args.append(arg)
            return new_args

        for i in merges_iter:
            if isinstance(instructions[n], startinput_class):
                instructions[n].args[1] += instructions[i].args[1]
            elif isinstance(instructions[n], (stopinput, gstopinput)):
                if instructions[n].get_size() != instructions[i].get_size():
                    raise NotImplemented()
                else:
                    instructions[n].args += instructions[i].args[1:]
            else:
                if instructions[n].get_size() != instructions[i].get_size():
                    # merge as non-vector instruction
                    instructions[n].args = expand_vector_args(instructions[n]) + \
                        expand_vector_args(instructions[i])
                    if instructions[n].is_vec():
                        instructions[n].size = 1
                else:
                    instructions[n].args += instructions[i].args
                
            # join arg_formats if not special iterators
            # if not isinstance(instructions[n].arg_format, (itertools.repeat, itertools.cycle)) and \
            #     not isinstance(instructions[i].arg_format, (itertools.repeat, itertools.cycle)):
            #     instructions[n].arg_format += instructions[i].arg_format
            instructions[i] = None
            self.merge_nodes(n, i)
            mergecount += 1

        return mergecount, n

    def compute_max_depths(self, depth_of):
        """ Compute the maximum 'depth' at which every instruction can be placed.
        This is the minimum depth of any merge_node succeeding an instruction.

        Similar to DAG shortest paths algorithm. Traverses the graph in reverse
        topological order, updating the max depth of each node's predecessors.
        """
        G = self.G
        merge_nodes_set = self.open_nodes
        top_order = Compiler.graph.topological_sort(G)
        max_depth_of = [None] * len(G)
        max_depth = max(depth_of)

        for i in range(len(max_depth_of)):
            if i in merge_nodes_set:
                max_depth_of[i] = depth_of[i] - 1
            else:
                max_depth_of[i] = max_depth

        for u in reversed(top_order):
            for v in G.pred[u]:
                if v not in merge_nodes_set:
                    max_depth_of[v] = min(max_depth_of[u], max_depth_of[v])
        return max_depth_of

    def merge_inputs(self):
        merges = defaultdict(list)
        remaining_input_nodes = []
        def do_merge(nodes):
            if len(nodes) > 1000:
                print 'Merging %d inputs...' % len(nodes)
            self.do_merge(iter(nodes))
        for n in self.input_nodes:
            inst = self.instructions[n]
            merge = merges[inst.args[0],inst.__class__]
            if len(merge) == 0:
                remaining_input_nodes.append(n)
            merge.append(n)
            if len(merge) >= self.max_parallel_open:
                do_merge(merge)
                merge[:] = []
        for merge in merges.itervalues():
            if merge:
                do_merge(merge)
        self.input_nodes = remaining_input_nodes

    def compute_preorder(self, merges, rev_depth_of):
        # find flexible nodes that can be on several levels
        # and find sources on level 0
        G = self.G
        merge_nodes_set = self.open_nodes
        depth_of = self.depths
        instructions = self.instructions
        flex_nodes = defaultdict(dict)
        starters = []
        for n in xrange(len(G)):
            if n not in merge_nodes_set and \
                depth_of[n] != rev_depth_of[n] and G[n] and G.get_attr(n,'start') == -1 and not isinstance(instructions[n], AsymmetricCommunicationInstruction):
                    #print n, depth_of[n], rev_depth_of[n]
                    flex_nodes[depth_of[n]].setdefault(rev_depth_of[n], set()).add(n)
            elif len(G.pred[n]) == 0 and \
                    not isinstance(self.instructions[n], RawInputInstruction):
                starters.append(n)
            if n % 10000000 == 0 and n > 0:
                print "Processed %d nodes at" % n, time.asctime()

        inputs = defaultdict(list)
        for node in self.input_nodes:
            player = self.instructions[node].args[0]
            inputs[player].append(node)
        first_inputs = [l[0] for l in inputs.itervalues()]
        other_inputs = []
        i = 0
        while True:
            i += 1
            found = False
            for l in inputs.itervalues():
                if i < len(l):
                    other_inputs.append(l[i])
                    found = True
            if not found:
                break
        other_inputs.reverse()

        preorder = []
        # magical preorder for topological search
        max_depth = max(merges)
        if max_depth > 10000:
            print "Computing pre-ordering ..."
        for i in xrange(max_depth, 0, -1):
            preorder.append(G.get_attr(merges[i], 'stop'))
            for j in flex_nodes[i-1].itervalues():
                preorder.extend(j)
            preorder.extend(flex_nodes[0].get(i, []))
            preorder.append(merges[i])
            if i % 100000 == 0 and i > 0:
                print "Done level %d at" % i, time.asctime()
        preorder.extend(other_inputs)
        preorder.extend(starters)
        preorder.extend(first_inputs)
        if max_depth > 10000:
            print "Done at", time.asctime()
        return preorder

    def compute_continuous_preorder(self, merges, rev_depth_of):
        print 'Computing pre-ordering for continuous computation...'
        preorder = []
        sources_for = defaultdict(list)
        stops_in = defaultdict(list)
        startinputs = []
        stopinputs = []
        for source in self.sources:
            sources_for[rev_depth_of[source]].append(source)
        for merge in merges.itervalues():
            stop = self.G.get_attr(merge, 'stop')
            stops_in[rev_depth_of[stop]].append(stop)
        for node in self.input_nodes:
            if isinstance(self.instructions[node], startinput_class):
                startinputs.append(node)
            else:
                stopinputs.append(node)
        max_round = max(rev_depth_of)
        for i in xrange(max_round, 0, -1):
            preorder.extend(reversed(stops_in[i]))
            preorder.extend(reversed(sources_for[i]))
        # inputs at the beginning
        preorder.extend(reversed(stopinputs))
        preorder.extend(reversed(sources_for[0]))
        preorder.extend(reversed(startinputs))
        return preorder

    def longest_paths_merge(self):
        """ Attempt to merge instructions of type instruction_type (which are given in
        merge_nodes) using longest paths algorithm.

        Returns the no. of rounds of communication required after merging (assuming 1 round/instruction).

        If reorder_between_opens is True, will attempt to place non-opens between start/stop opens.

        Doesn't use networkx.
        """
        G = self.G
        instructions = self.instructions
        #print(instructions)
        merge_nodes = self.open_nodes
        #print(merge_nodes)
        depths = self.depths
        #print(depths)

        if not merge_nodes and not self.input_nodes:
            return 0

        # merge opens at same depth
        merges = defaultdict(list)
        for node in merge_nodes:
            merges[depths[node]].append(node)

        # after merging, the first element in merges[i] remains for each depth i,
        # all others are removed from instructions and G
        last_nodes = [None, None]
        for i in sorted(merges):
            merge = merges[i]

            ### added for debug (start) ###
            # print("merge"+str(i)+":"+str(merge))
            ### added for debug (ended) ###

            if len(merge) > 1000:
                print 'Merging %d opens in round %d/%d' % (len(merge), i, len(merges))
            nodes = defaultdict(lambda: None)
            for b in (False, True):
                #my_merge = (m for m in merge if instructions[m] is not None and instructions[m].is_gf2n() is b)
                my_merge = (m for m in merge if instructions[m] is not None and isinstance(instructions[m], e_mult_class) is b)
                
                mc, nodes[0,b] = self.do_merge(iter(my_merge))

            # add edges to retain order of gf2n/modp start/stop opens
            for j in (0,1):
                node2 = nodes[j,True]
                nodep = nodes[j,False]
                if nodep is not None and node2 is not None:
                    G.add_edge(nodep, node2)
                # add edge to retain order of opens over rounds
                if last_nodes[j] is not None:
                    G.add_edge(last_nodes[j], node2 if nodep is None else nodep)
                last_nodes[j] = nodep if node2 is None else node2
            merges[i] = last_nodes[0]

        self.merge_inputs()

        preorder = None

        if len(instructions) > 100000:
            print "Topological sort ..."
        order = Compiler.graph.topological_sort(G, preorder)
        instructions[:] = [instructions[i] for i in order if instructions[i] is not None]
        if len(instructions) > 100000:
            print "Done at", time.asctime()

        return len(merges)

    def dependency_graph(self, merge_classes):
        """ Create the program dependency graph. """
        if len(merge_classes) != 1:
            if int(self.options.max_parallel_open):
                raise NotImplementedError('parallel limit only implemented ' \
                                          'for single instruction')

        block = self.block
        options = self.options
        open_nodes = set()
        self.open_nodes = open_nodes
        self.input_nodes = []
        colordict = defaultdict(lambda: 'gray', asm_open='red',\
                                ldi='lightblue', ldm='lightblue', stm='blue',\
                                mov='yellow', mulm='orange', mulc='orange',\
                                triple='green', square='green', bit='green',\
                                asm_input='lightgreen')

        G = Compiler.graph.SparseDiGraph(len(block.instructions))
        self.G = G

        reg_nodes = {}
        last_def = defaultdict(lambda: -1)
        last_mem_write = []
        last_mem_read = []
        warned_about_mem = []
        last_mem_write_of = defaultdict(list)
        last_mem_read_of = defaultdict(list)
        last_print_str = None
        last = defaultdict(lambda: defaultdict(lambda: None))
        last_open = deque()

        depths = [0] * len(block.instructions)
        self.depths = depths
        parallel_open = defaultdict(lambda: 0)
        next_available_depth = {}
        self.sources = []
        self.real_depths = [0] * len(block.instructions)
        round_type = {}

        def add_edge(i, j):
            G.add_edge(i, j)
            for d in (self.depths, self.real_depths):
                if d[j] < d[i]:
                    d[j] = d[i]

        def read(reg, n):
            if last_def[reg] != -1:
                add_edge(last_def[reg], n)

        def write(reg, n):
            last_def[reg] = n

        def handle_mem_access(addr, reg_type, last_access_this_kind,
                              last_access_other_kind):
            this = last_access_this_kind[addr,reg_type]
            other = last_access_other_kind[addr,reg_type]
            if this and other:
                if this[-1] < other[0]:
                    del this[:]
            this.append(n)
            for inst in other:
                add_edge(inst, n)

        def mem_access(n, instr, last_access_this_kind, last_access_other_kind):
            addr = instr.args[1]
            reg_type = instr.args[0].reg_type
            if isinstance(addr, int):
                for i in range(min(instr.get_size(), 100)):
                    addr_i = addr + i
                    handle_mem_access(addr_i, reg_type, last_access_this_kind,
                                      last_access_other_kind)
                if not warned_about_mem and (instr.get_size() > 100):
                    print 'WARNING: Order of memory instructions ' \
                        'not preserved due to long vector, errors possible'
                    warned_about_mem.append(True)
            else:
                handle_mem_access(addr, reg_type, last_access_this_kind,
                                  last_access_other_kind)
            if not warned_about_mem and not isinstance(instr, DirectMemoryInstruction):
                print 'WARNING: Order of memory instructions ' \
                    'not preserved, errors possible'
                # hack
                warned_about_mem.append(True)

        def keep_order(instr, n, t, arg_index=None):
            if arg_index is None:
                player = None
            else:
                player = instr.args[arg_index]
            if last[t][player] is not None:
                add_edge(last[t][player], n)
            last[t][player] = n

        for n,instr in enumerate(block.instructions):
            outputs,inputs = instr.get_def(), instr.get_used()

            G.add_node(n)

            # if options.debug:
            #     col = colordict[instr.__class__.__name__]
            #     G.add_node(n, color=col, label=str(instr))
            for reg in inputs:
                if reg.vector and instr.is_vec():
                    for i in reg.vector:
                        read(i, n)
                else:
                    read(reg, n)

            for reg in outputs:
                if reg.vector and instr.is_vec():
                    for i in reg.vector:
                        write(i, n)
                else:
                    write(reg, n)

            if isinstance(instr, merge_classes):
                open_nodes.add(n)
                G.add_node(n, merges=[])
                # the following must happen after adding the edge
                self.real_depths[n] += 1
                depth = depths[n] + 1
                while depth in round_type:
                    if round_type[depth] == type(instr):
                        break
                    depth += 1
                round_type[depth] = type(instr)
                if int(options.max_parallel_open):
                    skipped_depths = set()
                    while parallel_open[depth] >= int(options.max_parallel_open):
                        skipped_depths.add(depth)
                        depth = next_available_depth.get(depth, depth + 1)
                    for d in skipped_depths:
                        next_available_depth[d] = depth
                else:
                    self.real_depths[n] = depth
                parallel_open[depth] += len(instr.args) * instr.get_size()
                depths[n] = depth

            if isinstance(instr, ReadMemoryInstruction):
                if options.preserve_mem_order:
                    if last_mem_write and last_mem_read and last_mem_write[-1] > last_mem_read[-1]:
                        last_mem_read[:] = []
                    last_mem_read.append(n)
                    for i in last_mem_write:
                        add_edge(i, n)
                else:
                    mem_access(n, instr, last_mem_read_of, last_mem_write_of)
            elif isinstance(instr, WriteMemoryInstruction):
                if options.preserve_mem_order:
                    if last_mem_write and last_mem_read and last_mem_write[-1] < last_mem_read[-1]:
                        last_mem_write[:] = []
                    last_mem_write.append(n)
                    for i in last_mem_read:
                        add_edge(i, n)
                else:
                    mem_access(n, instr, last_mem_write_of, last_mem_read_of)
            # keep I/O instructions in order
            elif isinstance(instr, IOInstruction):
                if last_print_str is not None:
                    add_edge(last_print_str, n)
                last_print_str = n
            elif isinstance(instr, PublicFileIOInstruction):
                keep_order(instr, n, instr.__class__)
            elif isinstance(instr, RawInputInstruction):
                keep_order(instr, n, instr.__class__, 0)
                self.input_nodes.append(n)
                G.add_node(n, merges=[])
                player = instr.args[0]
                if isinstance(instr, stopinput):
                    add_edge(last[startinput_class][player], n)
                elif isinstance(instr, gstopinput):
                    add_edge(last[gstartinput][player], n)
            elif isinstance(instr, startprivateoutput_class):
                keep_order(instr, n, startprivateoutput_class, 2)
            elif isinstance(instr, stopprivateoutput_class):
                keep_order(instr, n, stopprivateoutput_class, 1)
            elif isinstance(instr, prep_class):
                keep_order(instr, n, instr.args[0])
            elif isinstance(instr, StackInstruction):
                keep_order(instr, n, StackInstruction)

            if not G.pred[n]:
                self.sources.append(n)

            if n % 100000 == 0 and n > 0:
                print "Processed dependency of %d/%d instructions at" % \
                    (n, len(block.instructions)), time.asctime()

        if len(open_nodes) > 1000:
            print "Program has %d %s instructions" % (len(open_nodes), merge_classes)

    def merge_nodes(self, i, j):
        """ Merge node j into i, removing node j """
        G = self.G
        if j in G[i]:
            G.remove_edge(i, j)
        if i in G[j]:
            G.remove_edge(j, i)
        G.add_edges_from(zip(itertools.cycle([i]), G[j], [G.weights[(j,k)] for k in G[j]]))
        G.add_edges_from(zip(G.pred[j], itertools.cycle([i]), [G.weights[(k,j)] for k in G.pred[j]]))
        G.get_attr(i, 'merges').append(j)
        G.remove_node(j)

    def eliminate_dead_code(self):
        instructions = self.instructions
        G = self.G
        merge_nodes = self.open_nodes
        count = 0
        open_count = 0
        for i,inst in zip(xrange(len(instructions) - 1, -1, -1), reversed(instructions)):
            # remove if instruction has result that isn't used
            unused_result = not G.degree(i) and len(inst.get_def()) \
                and reduce(operator.and_, (reg.can_eliminate for reg in inst.get_def())) \
                and not isinstance(inst, (DoNotEliminateInstruction))
            stop_node = G.get_attr(i, 'stop')
            unused_startopen = stop_node != -1 and instructions[stop_node] is None
            if unused_result or unused_startopen:
                G.remove_node(i)
                merge_nodes.discard(i)
                instructions[i] = None
                count += 1
                if unused_startopen:
                    open_count += len(inst.args)
        if count > 0:
            print 'Eliminated %d dead instructions, among which %d opens' % (count, open_count)

    def print_graph(self, filename):
        f = open(filename, 'w')
        print >>f, 'digraph G {'
        for i in range(self.G.n):
            for j in self.G[i]:
                print >>f, '"%d: %s" -> "%d: %s";' % \
                    (i, self.instructions[i], j, self.instructions[j])
        print >>f, '}'
        f.close()

    def print_depth(self, filename):
        f = open(filename, 'w')
        for i in range(self.G.n):
            print >>f, '%d: %s' % (self.depths[i], self.instructions[i])
        f.close()
