class Queue:
    def __init__(self, ):
        self._data = []

    def Que_in(self, e):
        self._data.append(e)

    def Que_out(self):
        return self._data.pop(0)

    def Que_first(self):
        return self._data[0]

    def Que_isEmpty(self):
        return len(self._data) == 0

    def Que_len(self):
        return len(self._data)

class ListNode:
    def __init__(self, value=0, next=None):
        self.value = value
        self.next = next


class LinkedList:
    def __init__(self):
        self.head = None
        self._length = 0  # 初始化长度计数器

    def __iter__(self):
        self._current = self.head
        return self

    def __next__(self):
        if self._current is None:
            raise StopIteration
        else:
            current_value = self._current.value
            self._current = self._current.next
            return current_value

    def append(self, value):
        if not self.head:
            self.head = ListNode(value)
        else:
            current = self.head
            while current.next:
                current = current.next
            current.next = ListNode(value)
        self._length += 1  # 链表长度加1

    def insert(self, value, position):
        if position < 0 or position > self._length:
            raise IndexError("Index out of bounds.")

        new_node = ListNode(value)
        if position == 0:
            new_node.next = self.head
            self.head = new_node
        else:
            current = self.head
            while position - 1:
                current = current.next
                position -= 1
            new_node.next = current.next
            current.next = new_node
        self._length += 1  # 链表长度加1

    def delete(self, value):
        current = self.head
        prev = None
        while current and current.value != value:
            prev = current
            current = current.next
        if current:  # 找到了要删除的节点
            if prev is None:
                self.head = current.next
            else:
                prev.next = current.next
            self._length -= 1  # 链表长度减1

    def find(self, value):
        current = self.head
        while current and current.value != value:
            current = current.next
        return current is not None

    def print_list(self):
        current = self.head
        while current:
            print(current.value, end=" -> ")
            current = current.next
        print("None")

    def to_list(self):
        node_values = []
        current = self.head
        while current:
            node_values.append(current.value)
            current = current.next
        return node_values

