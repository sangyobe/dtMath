/*!
\file      dtCscAMD.tpp
\brief      dtMath, Cholesky decomposition(L*L^T form) for Sparse Matrix Class
\author     Muhammad Zahak Jamal, zahakj@gmail.com
\author     Who is next author?
\date      2024. 06. 10
\version    1.0.0
\warning    Do Not delete this comment for document history! This is minimal manners!
*/

#ifndef DTMATH_DTCSC_AMD_TPP_
#define DTMATH_DTCSC_AMD_TPP_

#include "dtCscAMD.h"

namespace dt
{
namespace Math
{

/** Initializes the  CscAMD object from a  CscMatrix object
 */

template <uint16_t m_row, uint16_t m_col, typename m_type>
inline CscAMD<m_row, m_col, m_type>::CscAMD()
{
    m_elemNum = 0;
    memset(m_elem, 0, sizeof(m_type) * m_row * m_col);
    memset(m_rowIdx, 0, sizeof(int) * m_row * m_col);
    memset(m_colPtr, 0, sizeof(int) * (m_col + 1));
}

/** Initializes the  CscAMD object from a  CscMatrix object

* Function Arguments:
*
* m   ->   CscMatrix object
*/

template <uint16_t m_row, uint16_t m_col, typename m_type>
inline CscAMD<m_row, m_col, m_type>::CscAMD(const CscMatrix<m_row, m_col, m_type> &m)
{
    AMDOrder(m);
}

/** Finds the permutation vectors P and Pinv from a  CscMatrix object

* Function Arguments:
*
* m         ->   CscMatrix object
*
* \return    ->  Zero on success; -1 or failure
*/

template <uint16_t m_row, uint16_t m_col, typename m_type>
inline int CscAMD<m_row, m_col, m_type>::AMDOrder(const CscMatrix<m_row, m_col, m_type> &m)
{
    if (m_row != m_col)
    {
        assert(false && "Input matrix is not a square matrix");
        return -1;
    }

    m_elemNum = m.m_elemNum;
    memcpy(m_elem, m.m_elem, sizeof(m_type) * m_elemNum);
    memcpy(m_rowIdx, m.m_rowIdx, sizeof(int) * m_elemNum);
    memcpy(m_colPtr, m.m_colPtr, sizeof(int) * (m_col + 1));

    if (m_col != m_row)
        return -1;

    Compute();

    returnPinv();

    flag = 1;
    return 0;
}

/** Clears the W flag

*
*/

template <uint16_t m_row, uint16_t m_col, typename m_type>
inline int CscAMD<m_row, m_col, m_type>::wclear(int mark, int lemax, int *w)
{
    int k;
    if (mark < 2 || (mark + lemax < 0))
    {
        for (k = 0; k < m_col; k++)
        {
            if (w[k] != 0)
            {
                w[k] = 1;
            }
        }
        mark = 2; /* at this point, w [0..n-1] < mark holds */
    }
    return mark;
}

/** Finds the permutation vectors P by fixing the graph

*
*/

template <uint16_t m_row, uint16_t m_col, typename m_type>
inline int CscAMD<m_row, m_col, m_type>::tdfs(int j, int k, int *head, int *next, int *post, int *stck)
{
    int i, p, top = 0;
    if (!head || !next || !post || !stck) return (-1); /* check inputs */
    stck[0] = j;                                       /* place j on the stack */
    while (top >= 0)                                   /* while (stack is not empty) */
    {
        p = stck[top]; /* p = top of stack */
        i = head[p];   /* i = youngest child of p */
        if (i == -1)
        {
            top--;         /* p has no unordered children left */
            post[k++] = p; /* node p is the kth postordered node */
        }
        else
        {
            head[p] = next[i]; /* remove i from children of p */
            stck[++top] = i;   /* start dfs on child node i */
        }
    }
    return k;
}

/** Computes the permutation vectors P

*
*/

template <uint16_t m_row, uint16_t m_col, typename m_type>
inline int CscAMD<m_row, m_col, m_type>::Compute()
{
    int *last, nv[m_col + 1], w[m_col + 1], len[m_col + 1], next[m_col + 1], head[m_col + 1], elen[m_col + 1];
    int hhead[m_col + 1], degree[m_col + 1];

    int i, j, k, p, q, dk, mark, nel = 0, dense, d, mindeg = 0, elenk, nvi, nvj, nvk, nzmax, k3,
                                 pk1, pk2, k1, e, pj, k2, ln, lemax = 0, pk, eln, wnvi, p1, p2, p3, p4, pn, h, dext, jlast,
                                 ok, elemNum = m_elemNum;

    int rowIdx_s[m_row * m_col];
    int colPtr_s[m_row * m_col];

    memcpy(rowIdx_s, m_rowIdx, sizeof(int) * m_elemNum);
    memcpy(colPtr_s, m_colPtr, sizeof(int) * (m_col + 1));

    dense = std::max(16.0, 10 * sqrt((double)m_col));
    dense = std::min(m_col - 2, dense);

    nzmax = elemNum + elemNum / 5 + 2 * m_col;

    for (k = 0; k < m_col; k++)
    {
        len[k] = colPtr_s[k + 1] - colPtr_s[k];
    }
    len[m_col] = 0;
    last = P;
    /* --- Initialize quotient graph ---------------------------------------- */
    for (i = 0; i <= m_col; i++)
    {
        head[i] = -1; /* degree list i is empty */
        last[i] = -1;
        next[i] = -1;
        hhead[i] = -1;      /* hash list i is empty */
        nv[i] = 1;          /* node i is just one node --> In Paper i = {i} */
        w[i] = 1;           /* node i is alive */
        elen[i] = 0;        /* Ek of node i is empty */
        degree[i] = len[i]; /* degree of node i --> In Paper d[i] = |A[i]| */
    }

    mark = wclear(0, 0, w);

    /* --- Initialize degree lists ------------------------------------------ */
    for (i = 0; i < m_col; i++)
    {
        bool has_diag = false;
        for (p = colPtr_s[i]; p < colPtr_s[i + 1]; ++p)
        {
            if (rowIdx_s[p] == i)
            {
                has_diag = true;
                break;
            }
        }
        d = degree[i];

        if (d == 1 && has_diag) /* if node i is empty */
        {
            elen[i] = -2;     /* element i is dead */
            nel++;            /* add to number of elements in graph */
            colPtr_s[i] = -1; /* i is a root of assembly tree */
            w[i] = 0;
        }

        else if (d > dense || !has_diag) /* if node i is dense or has no structural diagonal element */
        {
            nv[i] = 0;    /* absorb i into element n */
            elen[i] = -1; /* node i is dead */
            nel++;        /* add to number of elements in graph */
            colPtr_s[i] = flip(m_col);
            nv[m_col]++;
        }
        else
        {
            if (head[d] != -1)
            {
                last[head[d]] = i;
            }
            next[i] = head[d]; /* put node i in degree list d */
            head[d] = i;
        }
    }

    colPtr_s[m_col] = -1;
    elen[m_col] = -2;
    w[m_col] = 0;

    while (nel < m_col) /* while (selecting pivots) do */
    {
        /* --- Select node of minimum approximate degree -------------------- */
        for (k = -1; mindeg < m_col && (k = head[mindeg]) == -1; mindeg++) {}
        if (next[k] != -1) last[next[k]] = -1;
        head[mindeg] = next[k]; /* remove k from degree list */
        elenk = elen[k];        /* elenk = |Ek| */
        nvk = nv[k];            /* # of nodes k represents */
        nel += nvk;             /* nv[k] nodes of A eliminated */

        /* --- Garbage collection ------------------------------------------- */
        if (elenk > 0 && elemNum + mindeg >= nzmax)
        {
            for (j = 0; j < m_col; j++)
            {
                if ((p = colPtr_s[j]) >= 0) /* j is a live node or element */
                {
                    colPtr_s[j] = rowIdx_s[p]; /* save first entry of object */
                    rowIdx_s[p] = flip(j);     /* first entry is now amd_flip(j) */
                }
            }
            for (q = 0, p = 0; p < elemNum;) /* scan all of memory */
            {
                if ((j = flip(rowIdx_s[p++])) >= 0) /* found object j */
                {
                    rowIdx_s[q] = colPtr_s[j]; /* restore first entry of object */
                    colPtr_s[j] = q++;         /* new pointer to object j */
                    for (k3 = 0; k3 < len[j] - 1; k3++)
                        rowIdx_s[q++] = rowIdx_s[p++];
                }
            }
            elemNum = q; /* Ci[m_elemNum...nzmax-1] now free */
        }

        /* --- Construct new element ---------------------------------------- */
        dk = 0;
        nv[k] = -nvk; /* flag k as in Lk */
        p = colPtr_s[k];
        pk1 = (elenk == 0) ? p : elemNum; /* do in place if elen[k] == 0 */
        pk2 = pk1;
        for (k1 = 1; k1 <= elenk + 1; k1++)
        {
            if (k1 > elenk)
            {
                e = k;               /* search the nodes in k */
                pj = p;              /* list of nodes starts at Ci[pj]*/
                ln = len[k] - elenk; /* length of list of nodes in k */
            }
            else
            {
                e = rowIdx_s[p++]; /* search the nodes in e */
                pj = colPtr_s[e];
                ln = len[e]; /* length of list of nodes in e */
            }
            for (k2 = 1; k2 <= ln; k2++)
            {
                i = rowIdx_s[pj++];
                if ((nvi = nv[i]) <= 0) continue; /* node i dead, or seen */
                dk += nvi;                        /* degree[Lk] += size of node i */
                nv[i] = -nvi;                     /* negate nv[i] to denote i in Lk*/
                rowIdx_s[pk2++] = i;              /* place i in Lk */
                if (next[i] != -1) last[next[i]] = last[i];
                if (last[i] != -1) /* remove i from degree list */
                {
                    next[last[i]] = next[i];
                }
                else
                {
                    head[degree[i]] = next[i];
                }
            }
            if (e != k)
            {
                colPtr_s[e] = flip(k); /* absorb e into k */
                w[e] = 0;              /* e is now a dead element */
            }
        }
        if (elenk != 0) elemNum = pk2; /* Ci[cnz...nzmax] is free */
        degree[k] = dk;                /* external degree of k - |Lk\i| */
        colPtr_s[k] = pk1;             /* element k is in Ci[pk1..pk2-1] */
        len[k] = pk2 - pk1;
        elen[k] = -2; /* k is now an element */

        /* --- Find set differences ----------------------------------------- */
        mark = wclear(mark, lemax, w); /* clear w if necessary */
        for (pk = pk1; pk < pk2; pk++) /* scan 1: find |Le\Lk| */
        {
            i = rowIdx_s[pk];
            if ((eln = elen[i]) <= 0) continue; /* skip if elen[i] empty */
            nvi = -nv[i];                       /* nv[i] was negated */
            wnvi = mark - nvi;
            for (p = colPtr_s[i]; p <= colPtr_s[i] + eln - 1; p++) /* scan Ei */
            {
                e = rowIdx_s[p];
                if (w[e] >= mark)
                {
                    w[e] -= nvi; /* decrement |Le\Lk| */
                }
                else if (w[e] != 0) /* ensure e is a live element */
                {
                    w[e] = degree[e] + wnvi; /* 1st time e seen in scan 1 */
                }
            }
        }

        /* --- Degree update ------------------------------------------------ */
        for (pk = pk1; pk < pk2; pk++) /* scan2: degree update */
        {
            i = rowIdx_s[pk]; /* consider node i in Lk */
            p1 = colPtr_s[i];
            p2 = p1 + elen[i] - 1;
            pn = p1;
            for (h = 0, d = 0, p = p1; p <= p2; p++) /* scan Ei */
            {
                e = rowIdx_s[p];
                if (w[e] != 0) /* e is an unabsorbed element */
                {
                    dext = w[e] - mark; /* dext = |Le\Lk| */
                    if (dext > 0)
                    {
                        d += dext;          /* sum up the set differences */
                        rowIdx_s[pn++] = e; /* keep e in Ei */
                        h += e;             /* compute the hash of node i */
                    }
                    else
                    {
                        colPtr_s[e] = flip(k); /* aggressive absorb. e->k */
                        w[e] = 0;              /* e is a dead element */
                    }
                }
            }
            elen[i] = pn - p1 + 1; /* elen[i] = |Ei| */
            p3 = pn;
            p4 = p1 + len[i];
            for (p = p2 + 1; p < p4; p++) /* prune edges in Ai */
            {
                j = rowIdx_s[p];
                if ((nvj = nv[j]) <= 0) continue; /* node j dead or in Lk */
                d += nvj;                         /* degree(i) += |j| */
                rowIdx_s[pn++] = j;               /* place j in node list of i */
                h += j;                           /* compute hash for node i */
            }
            if (d == 0) /* check for mass elimination */
            {
                colPtr_s[i] = flip(k); /* absorb i into k */
                nvi = -nv[i];
                dk -= nvi;  /* |Lk| -= |i| */
                nvk += nvi; /* |k| += nv[i] */
                nel += nvi;
                nv[i] = 0;
                elen[i] = -1; /* node i is dead */
            }
            else
            {
                degree[i] = std::min(degree[i], d); /* update degree(i) */
                rowIdx_s[pn] = rowIdx_s[p3];        /* move first node to end */
                rowIdx_s[p3] = rowIdx_s[p1];        /* move 1st el. to end of Ei */
                rowIdx_s[p1] = k;                   /* add k as 1st element in of Ei */
                len[i] = pn - p1 + 1;               /* new len of adj. list of node i */
                h %= m_col;                         /* finalize hash of i */
                next[i] = hhead[h];                 /* place i in hash bucket */
                hhead[h] = i;
                last[i] = h; /* save hash of i in last[i] */
            }
        } /* scan2 is done */
        degree[k] = dk; /* finalize |Lk| */
        lemax = std::max(lemax, dk);
        mark = wclear(mark + lemax, lemax, w); /* clear w */

        /* --- Supernode detection ------------------------------------------ */
        for (pk = pk1; pk < pk2; pk++)
        {
            i = rowIdx_s[pk];
            if (nv[i] >= 0) continue; /* skip if i is dead */
            h = last[i];              /* scan hash bucket of node i */
            i = hhead[h];
            hhead[h] = -1; /* hash bucket will be empty */
            for (; i != -1 && next[i] != -1; i = next[i], mark++)
            {
                ln = len[i];
                eln = elen[i];
                for (p = colPtr_s[i] + 1; p <= colPtr_s[i] + ln - 1; p++)
                    w[rowIdx_s[p]] = mark;
                jlast = i;
                for (j = next[i]; j != -1;) /* compare i with all j */
                {
                    ok = (len[j] == ln) && (elen[j] == eln);
                    for (p = colPtr_s[j] + 1; ok && p <= colPtr_s[j] + ln - 1; p++)
                    {
                        if (w[rowIdx_s[p]] != mark) ok = 0; /* compare i and j*/
                    }
                    if (ok) /* i and j are identical */
                    {
                        colPtr_s[j] = flip(i); /* absorb j into i */
                        nv[i] += nv[j];
                        nv[j] = 0;
                        elen[j] = -1; /* node j is dead */
                        j = next[j];  /* delete j from hash bucket */
                        next[jlast] = j;
                    }
                    else
                    {
                        jlast = j; /* j and i are different */
                        j = next[j];
                    }
                }
            }
        }

        /* --- Finalize new element------------------------------------------ */
        for (p = pk1, pk = pk1; pk < pk2; pk++) /* finalize Lk */
        {
            i = rowIdx_s[pk];
            if ((nvi = -nv[i]) <= 0) continue; /* skip if i is dead */
            nv[i] = nvi;                       /* restore nv[i] */
            d = degree[i] + dk - nvi;          /* compute external degree(i) */
            d = std::min(d, m_col - nel - nvi);
            if (head[d] != -1) last[head[d]] = i;
            next[i] = head[d]; /* put i back in degree list */
            last[i] = -1;
            head[d] = i;
            mindeg = std::min(mindeg, d); /* find new minimum degree */
            degree[i] = d;
            rowIdx_s[p++] = i; /* place i in Lk */
        }
        nv[k] = nvk;                 /* # nodes absorbed into k */
        if ((len[k] = p - pk1) == 0) /* length of adj list of element k*/
        {
            colPtr_s[k] = -1; /* k is a root of the tree */
            w[k] = 0;         /* k is now a dead element */
        }
        if (elenk != 0) elemNum = p; /* free unused space in Lk */
    }

    /* --- Postordering ----------------------------------------------------- */

    for (i = 0; i < m_col; i++)
        colPtr_s[i] = flip(colPtr_s[i]); /* fix assembly tree */
    for (j = 0; j <= m_col; j++)
        head[j] = -1;
    for (j = m_col; j >= 0; j--) /* place unordered nodes in lists */
    {
        if (nv[j] > 0) continue;     /* skip if j is an element */
        next[j] = head[colPtr_s[j]]; /* place j in list of its parent */
        head[colPtr_s[j]] = j;
    }
    for (e = m_col; e >= 0; e--) /* place elements in lists */
    {
        if (nv[e] <= 0) continue; /* skip unless e is an element */
        if (colPtr_s[e] != -1)
        {
            next[e] = head[colPtr_s[e]]; /* place e in list of its parent */
            head[colPtr_s[e]] = e;
        }
    }
    for (k = 0, i = 0; i <= m_col; i++) /* postorder the assembly tree */
    {
        if (colPtr_s[i] == -1) k = tdfs(i, k, head, next, P, w);
    }

    return 0;
}

/** Evaluates the Pinv vector

*
*/

template <uint16_t m_row, uint16_t m_col, typename m_type>
inline void CscAMD<m_row, m_col, m_type>::returnPinv()
{
    for (int j = 0; j < m_col; j++)
    {
        Pinv[P[j]] = j;
    }
}

/** Returns the P vector in  Vector

*
* \return    ->  Vector P in  Vector format
*/

template <uint16_t m_row, uint16_t m_col, typename m_type>
inline Vector<m_col, int> CscAMD<m_row, m_col, m_type>::GetP()
{
    if (!flag)
    {
        assert(false && "AMD ordering not complete");
        return Vector<m_col, int>();
    }

    return Vector<m_col, int>(P);
}

/** Returns the Pinv vector in  Vector

*
* \return    ->  Vector Pinv in  Vector format
*/

template <uint16_t m_row, uint16_t m_col, typename m_type>
inline Vector<m_col, int> CscAMD<m_row, m_col, m_type>::GetPinv()
{

    if (!flag)
    {
        assert(false && "AMD ordering not complete");
        return Vector<m_col, int>();
    }

    return Vector<m_col, int>(Pinv);
}

/** Permutes the input vector such that x = P * b

* Function Arguments:
*
* m         ->   CscMatrix object
*
* \return    ->  x = P * b in  Vector format
*/

template <uint16_t m_row, uint16_t m_col, typename m_type>
inline Vector<m_col, m_type> CscAMD<m_row, m_col, m_type>::PermuteVector(const Vector<m_col, m_type> &b)
{
    if (!flag)
    {
        assert(false && "AMD ordering not complete");
        return Vector<m_row, m_type>();
    }
    Vector<m_col, m_type> x;

    for (int k = 0; k < m_col; k++)
    {
        x.m_elem[Pinv[k]] = b.m_elem[k];
    }
    return Vector<m_row, m_type>(x);
}

/** Permutes the input vector such that x = P' * b

* Function Arguments:
*
* m         ->   CscMatrix object
*
* \return    ->  x = P' * b in  Vector format
*/

template <uint16_t m_row, uint16_t m_col, typename m_type>
inline Vector<m_col, m_type> CscAMD<m_row, m_col, m_type>::InvPermuteVector(const Vector<m_col, m_type> &b)
{
    if (!flag)
    {
        assert(false && "AMD ordering not complete");
        return Vector<m_col, m_type>();
    }

    Vector<m_col, m_type> x;

    for (int k = 0; k < m_col; k++)
    {
        x.m_elem[k] = b.m_elem[Pinv[k]];
    }
    return Vector<m_row, m_type>(x);
}

/** Permutes the matrix such that B = P * A * P'

*
* \return    ->  B = P * A * P' in  CscMatrix format
*/

template <uint16_t m_row, uint16_t m_col, typename m_type>
inline CscMatrix<m_row, m_col, m_type> CscAMD<m_row, m_col, m_type>::GetPermutedMatrix()
{

    if (!flag)
    {
        assert(false && "AMD ordering not complete");
        return CscMatrix<m_col, m_row, m_type>();
    }

    int cPtr, i, j, cnt = 0;
    m_type elem[m_elemNum];
    int colPtr[m_elemNum];
    int rowIdxCol[m_elemNum];

    for (j = 0; j < m_col; j++)
    {
        colPtr[j] = cnt;
        cPtr = P[j];
        for (i = m_colPtr[cPtr]; i < m_colPtr[cPtr + 1]; i++)
        {
            rowIdxCol[cnt] = Pinv[m_rowIdx[i]];
            elem[cnt++] = m_elem[i];
        }
    }

    colPtr[m_col] = cnt;

    return CscMatrix<m_col, m_row, m_type>(elem, m_elemNum, rowIdxCol, colPtr);
}

} // namespace Math
} // namespace dt

#endif // DTMATH_DTCSC_AMD_TPP_