
template <uint16_t mrow, uint16_t ncol, uint16_t order, typename m_type>
m_type CdjBezier<mrow, ncol, order, m_type>::BinomialCoeff(int n, int k)
{
    if (n < k) return 0;

    int min, max;
    m_type bc = 1;

    if ((n - k) >= k)
    {
        max = n - k;
        min = k;
    }
    else
    {
        max = k;
        min = n - k;
    }

    for (int i = 1; i <= min; i++)
    {
        bc *= ((m_type)(max + i) / (m_type)i);
    }

    return bc;
}

template <uint16_t mrow, uint16_t ncol, uint16_t order, typename m_type>
m_type CdjBezier<mrow, ncol, order, m_type>::BernsteinPoly(int n, int k, m_type ctrlParam)
{
    return BinomialCoeff(n, k) * std::pow(ctrlParam, k) * std::pow(1 - ctrlParam, n - k);
}

template <uint16_t mrow, uint16_t ncol, uint16_t order, typename m_type>
void CdjBezier<mrow, ncol, order, m_type>::BezierPoly(CdtMatrix<order, ncol, m_type> &ctrlPoints, m_type ctrlParam, m_type period)
{
    CdtMatrix<order - 1, ncol> mat1 = ctrlPoints.template GetBlock<order - 1, ncol>(0, 0);
    CdtMatrix<order - 1, ncol> mat2 = ctrlPoints.template GetBlock<order - 1, ncol>(1, 0);
    CdtMatrix<order - 1, ncol> delta_mat = mat2 - mat1;

    int n = order - 1;

    if (ctrlParam > 1)
    {
        ctrlParam = 1;
    }
    else if (ctrlParam < 0)
    {
        ctrlParam = 0;
    }

    pos.SetZero();
    vel.SetZero();
    for (int i = 0; i <= n; i++)
    {
        //pos_sum += BernsteinPoly(n, i, ctrlParam) * ctrlPoints.GetRowVec(i);
        pos += BernsteinPoly(n, i, ctrlParam) * ctrlPoints.GetRowVec(i);
        if (i <= n - 1)
        {
            //vel_sum += BernsteinPoly(n - 1, i, ctrlParam) * delta_mat.GetRowVec(i);
            vel += BernsteinPoly(n - 1, i, ctrlParam) * delta_mat.GetRowVec(i);
        }
    }

    //BezierOut<ncol> result;
    //result.pos = pos_sum;
    //if (period == 0) { period = 1; }
    //if (std::abs(period) <= std::numeric_limits<m_type>::epsilon()) period = 1;
    //result.vel = vel_sum * (m_type)n / period;

    //return result;

    if (std::abs(period) <= std::numeric_limits<m_type>::epsilon()) period = 1;
    vel = vel * (m_type)n / period;
}

template <uint16_t mrow, uint16_t ncol, uint16_t order, typename m_type>
void CdjBezier<mrow, ncol, order, m_type>::BezierInterp(CdtMatrix<mrow, ncol, m_type> &ctrlPoints, m_type ctrlParam, m_type period)
{
    CdtMatrix<mrow - 1, ncol> gapmat;
    m_type sum = 0;
    CdtVector<mrow> vec_index;

    gapmat = ctrlPoints.template GetBlock<mrow - 1, ncol>(1, 0) - ctrlPoints.template GetBlock<mrow - 1, ncol>(0, 0);
    vec_index.SetZero();

    for (int i = 0; i < mrow - 1; i++)
    {
        m_type temp_norm = gapmat.GetRowVec(i).GetNorm();
        vec_index(i + 1) = temp_norm + vec_index(i);
        sum += temp_norm;
    }

    //if (sum == 0) { sum = 1; }
    if (std::abs(sum) < std::numeric_limits<m_type>::epsilon()) sum = 1;
    vec_index /= sum;
    int index = 0;

    for (int i = 1; i < mrow - 1; i++)
    {
        if (ctrlParam >= vec_index(i)) { index++; }
        else { break; }
    }

    m_type ctrlPram_new = (ctrlParam - vec_index(index)) / (vec_index(index + 1) - vec_index(index)); // assigned control parameter
    CdtVector<mrow> period_set;
    period_set = vec_index * period;
    m_type period_new = (period_set(index + 1) - period_set(index)); // assigned period
    // index, ctrlPoints 
    //period_set.Print();
    //cout << index << endl;
    CdtMatrix<order, ncol> ctrlPoints_new;
    ctrlPoints_new.SetZero();

    CdtVector<order> designed_col;
    CdtVector<order / 2> vec1;
    CdtVector<order / 2> vec2;

    for (int i = 0; i < ncol; i++)
    {
        CdtVector<mrow> colvec = ctrlPoints.GetColVec(i);

        if (mrow <= 2)
        {
            vec1.SetFill(colvec(index));
            vec2.SetFill(colvec(index + 1));
            designed_col.SetBlock(0, vec1);
            designed_col.SetBlock(4, vec2);
        }
        else
        {
            if (index == 0)
            {
                vec1.SetFill(colvec(index));
                designed_col.SetBlock(0, vec1);

                if ((colvec(index + 1) - colvec(index)) * (colvec(index + 2) - colvec(index + 1)) > 0)
                {
                    m_type vecgap = colvec(index + 1) - colvec(index);
                    m_type period_next = (period_set(index + 2) - period_set(index + 1));
                    m_type vecgap_next = colvec(index + 2) - colvec(index + 1);

                    m_type delta_end = period_new / 14 * (vecgap / period_new + vecgap_next / period_next);

                    designed_col(4) = colvec(index + 1) - 3 * delta_end;
                    designed_col(5) = colvec(index + 1) - 2 * delta_end;
                    designed_col(6) = colvec(index + 1) - delta_end;
                    designed_col(7) = colvec(index + 1);

                }
                else
                {
                    vec2.SetFill(colvec(index + 1));
                    designed_col.SetBlock(4, vec2);
                }
            }
            else if (index == (mrow - 2))
            {
                if ((colvec(index + 1) - colvec(index)) * (colvec(index) - colvec(index - 1)) > 0)
                {
                    m_type vecgap = colvec(index + 1) - colvec(index);
                    m_type period_prev = (period_set(index) - period_set(index - 1));
                    m_type vecgap_prev = colvec(index) - colvec(index - 1);

                    m_type delta_start = period_new / 14 * (vecgap / period_new + vecgap_prev / period_prev);

                    designed_col(0) = colvec(index);
                    designed_col(1) = colvec(index) + delta_start;
                    designed_col(2) = colvec(index) + 2 * delta_start;
                    designed_col(3) = colvec(index) + 3 * delta_start;
                }
                else
                {
                    vec1.SetFill(colvec(index));
                    designed_col.SetBlock(0, vec1);
                }
                vec2.SetFill(colvec(index + 1));
                designed_col.SetBlock(4, vec2);
            }
            else
            {
                if ((colvec(index + 1) - colvec(index)) * (colvec(index) - colvec(index - 1)) > 0)
                {
                    m_type vecgap = colvec(index + 1) - colvec(index);
                    m_type period_prev = (period_set(index) - period_set(index - 1));
                    m_type vecgap_prev = colvec(index) - colvec(index - 1);

                    m_type delta_start = period_new / 14 * (vecgap / period_new + vecgap_prev / period_prev);

                    designed_col(0) = colvec(index);
                    designed_col(1) = colvec(index) + delta_start;
                    designed_col(2) = colvec(index) + 2 * delta_start;
                    designed_col(3) = colvec(index) + 3 * delta_start;
                }
                else
                {
                    vec1.SetFill(colvec(index));
                    designed_col.SetBlock(0, vec1);
                }
                vec2.SetFill(colvec(index + 1));
                designed_col.SetBlock(4, vec2);

                if ((colvec(index + 1) - colvec(index)) * (colvec(index + 2) - colvec(index + 1)) > 0)
                {
                    m_type vecgap = colvec(index + 1) - colvec(index);
                    m_type period_next = (period_set(index + 2) - period_set(index + 1));
                    m_type vecgap_next = colvec(index + 2) - colvec(index + 1);

                    m_type delta_end = period_new / 14 * (vecgap / period_new + vecgap_next / period_next);

                    designed_col(4) = colvec(index + 1) - 3 * delta_end;
                    designed_col(5) = colvec(index + 1) - 2 * delta_end;
                    designed_col(6) = colvec(index + 1) - delta_end;
                    designed_col(7) = colvec(index + 1);

                }
                else
                {
                    vec2.SetFill(colvec(index + 1));
                    designed_col.SetBlock(4, vec2);
                }
            }

        }
        //ctrlPoints_new.SetBlock<8, 1>(0, i, designed_col.GetElementsAddr());
        ctrlPoints_new.SetColVec(i, designed_col);
    }

    //Bezier<8, ncol> handler;
    //return handler.BezierPoly(ctrlPoints_new, ctrlPram_new, period_new);

    BezierPoly(ctrlPoints_new, ctrlPram_new, period_new);
}
