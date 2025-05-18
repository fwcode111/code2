end = -1


class Node(object):
    def __init__(self, leaf: bool):
        self.start = None
        self.end = None
        self.leaf = leaf
        self.children = {}
        self.suffix_link = None

    def __getattribute__(self, name: str):
        if name == 'end':
            if self.leaf:
                return end
        return super(Node, self).__getattribute__(name)


class suffixTree(object):
    def __init__(self, strings: str):
        self.T = strings + '#'
        self.actNode = None
        self.actEdge = -1
        self.actLength = 0
        self.remainder = 0
        self.insideNode = None
        self.len = -1
        self.root = None
        self.build_tree()

    def _edge_length(self, node: Node):
        return node.end - node.start + 1

    def _walk_down(self, node: Node):
        length = self._edge_length(node)
        if self.actLength >= length:
            self.actEdge += length
            self.actLength -= length
            self.actNode = node
            return True
        return False

    def _gen_node(self, start: int, end=None, leaf=False):
        node = Node(leaf)
        node.suffix_link = self.root
        node.start = start
        node.end = end
        return node

    def _gen_trie(self, pos: int):
        global end
        end = pos
        self.remainder += 1
        self.insideNode = None

        while self.remainder > 0:
            if self.actLength == 0:
                self.actEdge = pos

            if self.actNode.children.get(self.T[self.actEdge]) is None:
                self.actNode.children[self.T[self.actEdge]] = self._gen_node(pos, leaf=True)
                if self.insideNode:
                    self.insideNode.suffix_link = self.actNode
                    self.insideNode = None

            else:
                nextNode = self.actNode.children.get(self.T[self.actEdge])
                if self._walk_down(nextNode):
                    continue
                if self.T[nextNode.start + self.actLength] == self.T[pos]:
                    if self.insideNode and (self.actNode != self.root):
                        self.insideNode.suffix_link = self.actNode
                        self.insideNode = None

                    self.actLength += 1
                    break
                splitEnd = nextNode.start + self.actLength - 1
                split_node = self._gen_node(nextNode.start, splitEnd)
                self.actNode.children[self.T[self.actEdge]] = split_node
                split_node.children[self.T[pos]] = self._gen_node(pos, leaf=True)
                nextNode.start += self.actLength
                split_node.children[self.T[nextNode.start]] = nextNode

                if self.insideNode:
                    self.insideNode.suffix_link = split_node
                self.insideNode = split_node

            self.remainder -= 1
            if (self.actNode == self.root) and (self.remainder > 0):
                self.actLength -= 1
                self.actEdge = pos - self.remainder + 1
            elif self.actNode != self.root:
                self.actNode = self.actNode.suffix_link

    def build_tree(self):
        self.len = len(self.T)
        rootEnd = -1
        self.root = self._gen_node(-1, rootEnd)
        self.actNode = self.root
        for i in range(self.len): self._gen_trie(i)
