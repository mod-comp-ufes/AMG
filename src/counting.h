#include <stdlib.h>

#ifndef count_h
#define count_h

typedef struct Nodec {
    int i, m;
    struct Nodec *next, *prev;
} Nodec;

typedef struct {
    int n;
    Nodec *head, *tail;
} listc;

typedef struct {
    listc **count;
    Nodec **nodes;
    int *m;
    int min_m;
    int n;
} counting;

listc * create_lc() {
    listc *l = (listc *) malloc(sizeof (listc));
    l->n = 0;
    l->head = l->tail = NULL;
    return l;
};

void destroy_lc(listc *l) {
    Nodec *h;
    while (l->head) {
        h = l->head->next;
        free(l->head);
        l->head = h;
    }
    free(l);
};

counting * create_c(int n) {
    counting * c = (counting*) malloc(sizeof (counting));
    c->count = (listc**) malloc(2 * n * sizeof (listc*));
    for (int i = 0; i < 2 * n; i++) c->count[i] = create_lc();
    c->nodes = (Nodec**) malloc(n * sizeof (Nodec*));
    for (int i = 0; i < n; i++) c->nodes[i] = NULL;
    c->m = (int*) calloc(n, sizeof (int));
    c->n = n;
    c->min_m = 2 * n;
    return c;
}

void destroy_c(counting * c) {
    for (int i = 0; i < c->n; i++) destroy_lc(c->count[i]);
    free(c->count);
    free(c->nodes);
    free(c->m);
}

Nodec * insert_c(counting *c, int i, int m) {
    listc *l = c->count[m];
    Nodec *t = (Nodec *) malloc(sizeof (Nodec));
    c->nodes[i] = t;
    t->i = i;
    t->m = m;
    t->next = NULL;
    t->prev = l->tail;
    if (l->tail) {
        l->tail->next = t;
        l->tail = t;
    } else {
        l->head = l->tail = t;
    }
    l->n++;
    c->m[i] = m;
    if (c->m[i] < c->min_m)
        c->min_m = c->m[i];

    return t;
};

void remove_c(counting *c, int i) {
    Nodec *n = c->nodes[i];
    int m = n->m;

    if (n->prev) {
        n->prev->next = n->next;
    } else {
        c->count[m]->head = n->next;
    }

    if (n->next) {
        n->next->prev = n->prev;
    } else {
        c->count[m]->tail = n->prev;
    }

    free(n);
    c->nodes[i] = NULL;
    c->count[m]->n--;

    if (c->count[c->min_m]->n == 0) {
        for (c->min_m = m + 1; ((c->min_m < 2*(c->n)) && (c->count[c->min_m]->n == 0)); c->min_m++) {
        }
    }
}

Nodec* move_c(counting *c, int i, int new_m) {
    /* declarations */
    listc *l = c->count[new_m]; // new list
    Nodec *t = c->nodes[i]; // moving node
    int m = t->m; // old t

    // update m[] array
    c->m[i] = new_m;

    // remove t from old array
    if (t->prev) {
        t->prev->next = t->next;
    } else {
        c->count[m]->head = t->next;
    }
    if (t->next) {
        t->next->prev = t->prev;
    } else {
        c->count[m]->tail = t->prev;
    }

    // add t to new array
    t->m = new_m;
    t->next = NULL;
    t->prev = l->tail;
    if (l->tail) {
        l->tail->next = t;
        l->tail = t;
    } else {
        l->head = l->tail = t;
    }
    l->n++;

    // update nodes[] array
    c->nodes[i] = t;
    c->count[m]->n--;

    // update min_m
    if (new_m < c->min_m)
        c->min_m = new_m;
    else {
        if (c->count[c->min_m]->n == 0) {
            for (c->min_m = m + 1; ((c->min_m < 2*(c->n)) && (c->count[c->min_m]->n == 0)); c->min_m++) {
            }
        }
    }

    return t;
}

#endif