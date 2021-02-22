# include <stdlib.h>
# include "list.h"
list * create_l () {
	list *l = (list *) malloc(sizeof(list));
	l->n = 0;
	l->head = l->tail = NULL;
	return l;
};

void destroy_l (list *l) {
	Node *h;
	while (l->head) {
		h = l->head->next;
		free(l->head);
		l->head = h;
	}
	free(l);
};

void insert_l_tail (list *l, int e, double v) {
	Node *t = (Node *) malloc(sizeof(Node));
	t->elem = e; t->val = v; t->next = NULL;
	if (l->tail) {
		l->tail->next = t;
		l->tail = t;
	} else {
		l->head = l->tail = t;
	}
	l->n++;
};

void insert_l_head (list *l, int e, double v) {
	Node *t = (Node *) malloc(sizeof(Node));
	t->elem = e; t->val = v; 
	t->next = l->head;
	l->head = t;
	l->n++;
};

void insert_l_sort (list *l, int e, double v) { // ascending order of elem
	Node *t = (Node *) malloc(sizeof(Node)), *curr;
	t->elem = e; t->val = v;
	if (!l->n || e < l->head->elem) {
		t->next = l->head;
		l->head = t; if (!l->tail) l->tail = t;
	} else {
		for (curr=l->head; curr->next && e >= curr->next->elem; curr=curr->next) {}
		t->next = curr->next;
		curr->next = t;
		if (!t->next) l->tail = t;
	}
	l->n++;
};

Node* check_in_l(list *l, int i) {
	Node* aux;
	for(aux=l->head;((aux) && (aux->val != i));aux=aux->next) {}
	return aux;
}

list1 * create_l1 () {
	list1 *l = (list1 *) malloc(sizeof(list1));
	l->n = 0;
	l->head = NULL;
	return l;
};

void destroy_l1 (list1 *l) {
	Node1 *h;
	while (l->head) {
		h = l->head->next;
		free(l->head);
		l->head = h;
	}
	free(l);
};

void insert_l1 (list1 *l, int e) {
	Node1 *h = (Node1 *) malloc(sizeof(Node1));
	h->elem = e; h->next = l->head;
	l->head = h;
	l->n++;
};

int remove_fst (list1 *l) {
	int fst = l->head->elem;
	Node1 *h = l->head->next;
	free(l->head); l->head = h; l->n--;
	return fst;
};

list_p * create_l_p () {
	list_p *l = (list_p *) malloc(sizeof(list_p));
	l->n = 0;
	l->head = l->tail = NULL;
	return l;
};

void destroy_l_p (list_p* l) {
	Node_p *h;
	while (l->head) {
		h = l->head->next;
		free(l->head);
		l->head = h;
	}
	free(l);
};

void insert_l_p_tail (list_p* l, void* val, double elem) {
	Node_p *t = (Node_p *) malloc(sizeof(Node_p));
	t->val= val; t->elem=elem; t->next=NULL;
	if (l->tail) {
		l->tail->next = t;
		l->tail = t;
	} else {
		l->head = l->tail = t;
	}
	l->n++;
};
