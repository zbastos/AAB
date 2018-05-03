
class Trie:
    def __init__(self):
        self.nodes = {0: {}}  # root node
        self.num = 0

    def print_tree(self):
        for k in self.nodes.keys():
            print(k, "->", self.nodes[k])

    def addNode(self, origin, symbol):
        self.num += 1
        self.nodes[origin][symbol] = self.num
        self.nodes[self.num] = {}

    def addPattern(self, p):
        pos = 0
        node = 0
        while pos < len(p):
            if p[pos] not in self.nodes[node].keys():
                self.addNode(node, p[pos])
            node = self.nodes[node][p[pos]]
            pos += 1

    def trieFromPatterns(self, pats):
        for p in pats:
            self.addPattern(p)

    def prefixTrieMatch(self, text):
        pos = 0
        match = ""
        node = 0
        while pos < len(text):
            if text[pos] in self.nodes[node].keys():
                node = self.nodes[node][text[pos]]
                match += text[pos]
                if self.nodes[node] == {}:
                    return match
                else:
                    pos += 1
            else:
                return None
        return None

    def trieMatches(self, text):
        res = []
        for i in range(len(text)):
            m = self.prefixTrieMatch(text[i:])
            if m != None: res.append((i, m))
        return res


class SuffixTrie:

    def __init__(self):
        self.nodes = {0: (-1, {})}  # root node
        self.num = 0

    def print_tree(self):
        for k in self.nodes.keys():
            if self.nodes[k][0] < 0:
                print(k, "->", self.nodes[k][1])
            else:
                print(k, ":", self.nodes[k][0])

    def addNode(self, origin, symbol, leafnum=-1):
        self.num += 1
        self.nodes[origin][1][symbol] = self.num
        self.nodes[self.num] = (leafnum, {})

    def addSuffix(self, p, sufnum):
        pos = 0
        node = 0
        while pos < len(p):
            if p[pos] not in self.nodes[node][1].keys():
                if pos == len(p) - 1:
                    self.addNode(node, p[pos], sufnum)
                else:
                    self.addNode(node, p[pos])
            node = self.nodes[node][1][p[pos]]
            pos += 1

    def suffixTrieFromSeq(self, text):
        t = text + "$"
        for i in range(len(t)):
            self.addSuffix(t[i:], i)

    def findPattern(self, pattern):
        pos = 0
        node = 0
        for pos in range(len(pattern)):
            if pattern[pos] in self.nodes[node][1].keys():
                node = self.nodes[node][1][pattern[pos]]
                pos += 1
            else:
                return None
        return self.getLeafesBelow(node)

    def getLeafesBelow(self, node):
        res = []
        if self.nodes[node][0] >= 0:
            res.append(self.nodes[node][0])
        else:
            for k in self.nodes[node][1].keys():
                newnode = self.nodes[node][1][k]
                leafes = self.getLeafesBelow(newnode)
                res.extend(leafes)
        return res


if __name__=="__main__":
    patterns = ["ATAGA", "ATC", "GAT"]
    t = Trie()
    t.trieFromPatterns(patterns)
    t.print_tree()
    print (t.prefixTrieMatch("ATAGACATC"))
    print (t.trieMatches("CCATAGACATCAAGATCGG"))

    print ("-------- Suffix Trie ------")
    seq = "TACTA"
    st = SuffixTrie()
    st.suffixTrieFromSeq(seq)
    st.print_tree()
    print(st.findPattern("TA"))

