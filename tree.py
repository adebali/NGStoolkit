from ete3 import Tree

t = Tree('(t);')
a = t.add_child(Tree('(a);'))



a.add_sister(Tree('(s);'))
t.add_child(Tree('(ts);'))
b = a.add_child(Tree('(b);'))
c = b.add_child(Tree('(c);'))
d = b.add_child(Tree('(d);'))

print(b)

# class Node(object):
#     def __init__(self, input = ''):
#         self.parent = None
#         self.children = []
#         self.name = input

#     def get_children(self):
#         return self.children
    
#     def get_parent(self):
#         return self.parent

#     def get_name(self):
#         return self.name

#     def add_child(self, input = ''):
#         child = Node(input)
#         child.parent = self
#         self.children.append(child)


# mytree = Tree('firstNode')
# mytree.add_child('child1')
# print(mytree.get_children()[0].get_parent().get_name())