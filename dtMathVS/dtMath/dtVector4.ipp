#pragma once

template <typename m_type, uint16_t m_row>
inline CdtVector4<m_type, m_row>::CdtVector4()
{
    m_elem[0] = 0;
    m_elem[1] = 0;
    m_elem[2] = 0;
    m_elem[3] = 0;
}

template <typename m_type, uint16_t m_row>
inline CdtVector4<m_type, m_row>::CdtVector4(const m_type* element)
{
    m_elem[0] = element[0];
    m_elem[1] = element[1];
    m_elem[2] = element[2];
    m_elem[3] = element[3];
}

template <typename m_type, uint16_t m_row>
inline CdtVector4<m_type, m_row>::CdtVector4(const m_type* element, const size_t n_byte)
{
    switch (n_byte / sizeof(m_type))
    {
    case 1:
        m_elem[0] = element[0];
        m_elem[1] = 0;
        m_elem[2] = 0;
        m_elem[3] = 0;
        break;
    case 2:
        m_elem[0] = element[0];
        m_elem[1] = element[1];
        m_elem[2] = 0;
        m_elem[3] = 0;
        break;
    case 3:
        m_elem[0] = element[0];
        m_elem[1] = element[1];
        m_elem[2] = element[2];
        m_elem[3] = 0;
        break;
    default:
        m_elem[0] = element[0];
        m_elem[1] = element[1];
        m_elem[2] = element[2];
        m_elem[3] = element[3];
        break;
    }
}

template <typename m_type, uint16_t m_row>
inline CdtVector4<m_type, m_row>::CdtVector4(const m_type v0, const m_type v1, const m_type v2, const m_type v3)
{
    m_elem[0] = v0;
    m_elem[1] = v1;
    m_elem[2] = v2;
    m_elem[3] = v3;
}

template <typename m_type, uint16_t m_row>
inline CdtVector4<m_type, m_row>::CdtVector4(const CdtVector4& v)
{
    m_elem[0] = v.m_elem[0];
    m_elem[1] = v.m_elem[1];
    m_elem[2] = v.m_elem[2];
    m_elem[3] = v.m_elem[3];
}

template <typename m_type, uint16_t m_row>
inline CdtVector4<m_type, m_row>::CdtVector4(const CdtVector<m_row, m_type>& v)
{
    m_elem[0] = v.m_elem[0];
    m_elem[1] = v.m_elem[1];
    m_elem[2] = v.m_elem[2];
    m_elem[3] = v.m_elem[3];
}

template<typename m_type, uint16_t m_row>
inline CdtVector4<m_type, m_row>::CdtVector4(const CdtMatrix<m_row, 1, m_type>& v)
{
    m_elem[0] = v.m_elem[0];
    m_elem[1] = v.m_elem[1];
    m_elem[2] = v.m_elem[2];
    m_elem[3] = v.m_elem[3];
}

template <typename m_type, uint16_t m_row>
inline void CdtVector4<m_type, m_row>::SetZero()
{
    m_elem[0] = 0;
    m_elem[1] = 0;
    m_elem[2] = 0;
    m_elem[3] = 0;
}

template <typename m_type, uint16_t m_row>
inline void CdtVector4<m_type, m_row>::SetFill(const m_type value)
{
    m_elem[0] = value;
    m_elem[1] = value;
    m_elem[2] = value;
    m_elem[3] = value;
}

template <typename m_type, uint16_t m_row>
inline void CdtVector4<m_type, m_row>::SetElement(const m_type* element, const size_t n_byte)
{
    switch (n_byte / sizeof(m_type))
    {
    case 1:
        m_elem[0] = element[0];
        break;
    case 2:
        m_elem[0] = element[0];
        m_elem[1] = element[1];
        break;
    case 3:
        m_elem[0] = element[0];
        m_elem[1] = element[1];
        m_elem[2] = element[2];
        break;
    default:
        m_elem[0] = element[0];
        m_elem[1] = element[1];
        m_elem[2] = element[2];
        m_elem[3] = element[3];
        break;
    }
}

template <typename m_type, uint16_t m_row>
inline void CdtVector4<m_type, m_row>::SetElement(const m_type v0, const m_type v1, const m_type v2, const m_type v3)
{
    m_elem[0] = v0;
    m_elem[1] = v1;
    m_elem[2] = v2;
    m_elem[3] = v3;
}

template <typename m_type, uint16_t m_row>
inline void CdtVector4<m_type, m_row>::SetElement(const CdtVector4& v)
{
    m_elem[0] = v.m_elem[0];
    m_elem[1] = v.m_elem[1];
    m_elem[2] = v.m_elem[2];
    m_elem[3] = v.m_elem[3];
}

template <typename m_type, uint16_t m_row>
inline void CdtVector4<m_type, m_row>::SetElement(const CdtVector<m_row, m_type>& v)
{
    m_elem[0] = v.m_elem[0];
    m_elem[1] = v.m_elem[1];
    m_elem[2] = v.m_elem[2];
    m_elem[3] = v.m_elem[3];
}

template<typename m_type, uint16_t m_row>
inline void CdtVector4<m_type, m_row>::SetElement(const CdtMatrix<m_row, 1, m_type>& v)
{
    m_elem[0] = v.m_elem[0];
    m_elem[1] = v.m_elem[1];
    m_elem[2] = v.m_elem[2];
    m_elem[3] = v.m_elem[3];
}

template<typename m_type, uint16_t m_row>
template<uint16_t row>
inline void CdtVector4<m_type, m_row>::SetBlock(const uint16_t idxRow, const CdtVector<row, m_type>& v)
{
    if (idxRow >= m_row) return;

    uint16_t rowSz = m_row - idxRow;
    if (rowSz > row) rowSz = row;

    switch (rowSz)
    {
    case 1:
        m_elem[idxRow] = v.m_elem[0];
        break;
    case 2:
        m_elem[idxRow] = v.m_elem[0];
        m_elem[idxRow + 1] = v.m_elem[1];
        break;
    case 3:
        m_elem[idxRow] = v.m_elem[0];
        m_elem[idxRow + 1] = v.m_elem[1];
        m_elem[idxRow + 2] = v.m_elem[2];
        break;
    default:
        m_elem[0] = v.m_elem[0];
        m_elem[1] = v.m_elem[1];
        m_elem[2] = v.m_elem[2];
        m_elem[3] = v.m_elem[3];
        break;
    }
}

template<typename m_type, uint16_t m_row>
inline void CdtVector4<m_type, m_row>::SetBlock(const uint16_t idxRow, const m_type * v, const size_t n_byte)
{
    if (idxRow >= m_row) return;

    uint16_t rowSz = m_row - idxRow;
    uint16_t row = n_byte / sizeof(m_type);
    if (rowSz > row) rowSz = row;

    switch (rowSz)
    {
    case 1:
        m_elem[idxRow] = v[0];
        break;
    case 2:
        m_elem[idxRow] = v[0];
        m_elem[idxRow + 1] = v[1];
        break;
    case 3:
        m_elem[idxRow] = v[0];
        m_elem[idxRow + 1] = v[1];
        m_elem[idxRow + 2] = v[2];
        break;
    default:
        m_elem[0] = v[0];
        m_elem[1] = v[1];
        m_elem[2] = v[2];
        m_elem[3] = v[3];
        break;
    }
}

template<typename m_type, uint16_t m_row>
inline void CdtVector4<m_type, m_row>::SetBlock(const uint16_t idxRow, const CdtVector3<m_type, 3>& v)
{
    if (idxRow >= m_row) return;

    uint16_t rowSz = m_row - idxRow;
    if (rowSz > 3) rowSz = 3;

    switch (rowSz)
    {
    case 1:
        m_elem[idxRow] = v.m_elem[0];
        break;
    case 2:
        m_elem[idxRow] = v.m_elem[0];
        m_elem[idxRow + 1] = v.m_elem[1];
        break;
    case 3:
        m_elem[idxRow] = v.m_elem[0];
        m_elem[idxRow + 1] = v.m_elem[1];
        m_elem[idxRow + 2] = v.m_elem[2];
        break;
    }
}

template<typename m_type, uint16_t m_row>
inline void CdtVector4<m_type, m_row>::SetBlock(const uint16_t idxRow, const CdtVector4<m_type, 4>& v)
{
    if (idxRow >= m_row) return;

    switch (idxRow)
    {
    case 3:
        m_elem[3] = v.m_elem[0];
        break;
    case 2:
        m_elem[2] = v.m_elem[0];
        m_elem[3] = v.m_elem[1];
        break;
    case 1:
        m_elem[1] = v.m_elem[0];
        m_elem[2] = v.m_elem[1];
        m_elem[3] = v.m_elem[2];
        break;
    default:
        m_elem[0] = v.m_elem[0];
        m_elem[1] = v.m_elem[1];
        m_elem[2] = v.m_elem[2];
        m_elem[3] = v.m_elem[3];
        break;
    }
}

template<typename m_type, uint16_t m_row>
inline void CdtVector4<m_type, m_row>::SetBlock(const uint16_t idxRow, const CdtVector6<m_type, 6>& v)
{
    if (idxRow >= m_row) return;

    switch (idxRow)
    {
    case 3:
        m_elem[3] = v.m_elem[0];
        break;
    case 2:
        m_elem[2] = v.m_elem[0];
        m_elem[3] = v.m_elem[1];
        break;
    case 1:
        m_elem[1] = v.m_elem[0];
        m_elem[2] = v.m_elem[1];
        m_elem[3] = v.m_elem[2];
        break;
    default:
        m_elem[0] = v.m_elem[0];
        m_elem[1] = v.m_elem[1];
        m_elem[2] = v.m_elem[2];
        m_elem[3] = v.m_elem[3];
        break;
    }
}

template<typename m_type, uint16_t m_row>
template<uint16_t row>
inline void CdtVector4<m_type, m_row>::SetBlock(const uint16_t idxRow, const CdtMatrix<row, 1, m_type>& v)
{
    if (idxRow >= m_row) return;

    uint16_t rowSz = m_row - idxRow;
    if (rowSz > row) rowSz = row;

    switch (rowSz)
    {
    case 1:
        m_elem[idxRow] = v.m_elem[0];
        break;
    case 2:
        m_elem[idxRow] = v.m_elem[0];
        m_elem[idxRow + 1] = v.m_elem[1];
        break;
    case 3:
        m_elem[idxRow] = v.m_elem[0];
        m_elem[idxRow + 1] = v.m_elem[1];
        m_elem[idxRow + 2] = v.m_elem[2];
        break;
    default:
        m_elem[0] = v.m_elem[0];
        m_elem[1] = v.m_elem[1];
        m_elem[2] = v.m_elem[2];
        m_elem[3] = v.m_elem[3];
        break;
    }
}

template <typename m_type, uint16_t m_row>
inline void CdtVector4<m_type, m_row>::SetSwap(const uint16_t i, const uint16_t j)
{
    m_type elem = m_elem[i];
    m_elem[i] = m_elem[j];
    m_elem[j] = elem;
}

template <typename m_type, uint16_t m_row>
inline void CdtVector4<m_type, m_row>::SetNormalize()
{
    m_type norm = std::sqrt(
        m_elem[0] * m_elem[0] +
        m_elem[1] * m_elem[1] +
        m_elem[2] * m_elem[2] +
        m_elem[3] * m_elem[3]);

    if (norm < std::numeric_limits<m_type>::epsilon())
        norm = std::numeric_limits<m_type>::epsilon();

    m_elem[0] /= norm;
    m_elem[1] /= norm;
    m_elem[2] /= norm;
    m_elem[3] /= norm;
}

template <typename m_type, uint16_t m_row>
inline const m_type* const CdtVector4<m_type, m_row>::GetElementsAddr() const
{
    return m_elem;
}

template<typename m_type, uint16_t m_row>
template<uint16_t row>
inline CdtVector<row, m_type> CdtVector4<m_type, m_row>::GetBlock(const uint16_t idx)
{
    m_type elem[row] = { 0, };
    uint16_t rowSize = m_row - idx;

    if (idx >= m_row) return CdtVector<row, m_type>(elem);
    if (rowSize > row) rowSize = row;

    memcpy(elem, &m_elem[idx], sizeof(m_type) * rowSize);

    return CdtVector<row, m_type>(elem);
}

template<typename m_type, uint16_t m_row>
inline CdtVector3<m_type, 3> CdtVector4<m_type, m_row>::GetBlockVec3(const uint16_t idx)
{
    m_type elem[3] = { 0, };

    if (idx >= m_row) return CdtVector3<m_type, 3>(elem);

    switch (m_row - idx)
    {
    case 1:
        elem[0] = m_elem[idx];
        break;
    case 2:
        elem[0] = m_elem[idx];
        elem[1] = m_elem[idx + 1];
        break;
    default:
        elem[0] = m_elem[idx];
        elem[1] = m_elem[idx + 1];
        elem[2] = m_elem[idx + 2];
    };

    return CdtVector3<m_type, 3>(elem);
}

template<typename m_type, uint16_t m_row>
template<uint16_t row>
inline int8_t CdtVector4<m_type, m_row>::GetBlock(const uint16_t idx, CdtVector<row, m_type>& v)
{
    uint16_t rowSize = m_row - idx;

    if (idx >= m_row) return -1;
    if (rowSize > row) rowSize = row;

    memcpy(v.m_elem, &m_elem[idx], sizeof(m_type) * rowSize);

    return 0;
}

template<typename m_type, uint16_t m_row>
inline int8_t CdtVector4<m_type, m_row>::GetBlockVec3(const uint16_t idx, CdtVector3<m_type, 3>& v)
{
    if (idx >= m_row) return -1;

    switch (m_row - idx)
    {
    case 1:
        v.m_elem[0] = m_elem[idx];
        break;
    case 2:
        v.m_elem[0] = m_elem[idx];
        v.m_elem[1] = m_elem[idx + 1];
        break;
    default:
        v.m_elem[0] = m_elem[idx];
        v.m_elem[1] = m_elem[idx + 1];
        v.m_elem[2] = m_elem[idx + 2];
    };

    return 0;
}

template <typename m_type, uint16_t m_row>
inline m_type CdtVector4<m_type, m_row>::GetNorm() const
{
    return std::sqrt(
        m_elem[0] * m_elem[0] +
        m_elem[1] * m_elem[1] +
        m_elem[2] * m_elem[2] +
        m_elem[3] * m_elem[3]);
}

template <typename m_type, uint16_t m_row>
inline m_type CdtVector4<m_type, m_row>::GetSqNorm() const
{
    return (
        m_elem[0] * m_elem[0] +
        m_elem[1] * m_elem[1] +
        m_elem[2] * m_elem[2] +
        m_elem[3] * m_elem[3]);
}

template <typename m_type, uint16_t m_row>
inline m_type CdtVector4<m_type, m_row>::GetSum() const
{
    return (
        m_elem[0] +
        m_elem[1] +
        m_elem[2] +
        m_elem[3]);
}

template <typename m_type, uint16_t m_row>
inline CdtVector4<m_type, m_row> CdtVector4<m_type, m_row>::GetNormalized() const
{
    m_type norm = std::sqrt(
        m_elem[0] * m_elem[0] +
        m_elem[1] * m_elem[1] +
        m_elem[2] * m_elem[2] +
        m_elem[3] * m_elem[3]);

    if (norm < std::numeric_limits<m_type>::epsilon())
        norm = std::numeric_limits<m_type>::epsilon();

    return CdtVector4(
        m_elem[0] / norm,
        m_elem[1] / norm,
        m_elem[2] / norm,
        m_elem[3] / norm);
}

template <typename m_type, uint16_t m_row>
inline CdtMatrix<1, m_row, m_type> CdtVector4<m_type, m_row>::Transpose() const
{
    return CdtMatrix<1, m_row, m_type>(m_elem);
}

/* Assignment operators */
template <typename m_type, uint16_t m_row>
inline CdtVector4<m_type, m_row>& CdtVector4<m_type, m_row>::operator =(const CdtVector4& v)
{
    m_elem[0] = v.m_elem[0];
    m_elem[1] = v.m_elem[1];
    m_elem[2] = v.m_elem[2];
    m_elem[3] = v.m_elem[3];

    return (*this);
}

template <typename m_type, uint16_t m_row>
inline CdtVector4<m_type, m_row>& CdtVector4<m_type, m_row>::operator +=(const CdtVector4& v)
{
    m_elem[0] += v.m_elem[0];
    m_elem[1] += v.m_elem[1];
    m_elem[2] += v.m_elem[2];
    m_elem[3] += v.m_elem[3];

    return (*this);
}

template <typename m_type, uint16_t m_row>
inline CdtVector4<m_type, m_row>& CdtVector4<m_type, m_row>::operator -=(const CdtVector4& v)
{
    m_elem[0] -= v.m_elem[0];
    m_elem[1] -= v.m_elem[1];
    m_elem[2] -= v.m_elem[2];
    m_elem[3] -= v.m_elem[3];

    return (*this);
}

template <typename m_type, uint16_t m_row>
inline CdtVector4<m_type, m_row>& CdtVector4<m_type, m_row>::operator *=(const CdtVector4& v)
{
    m_elem[0] *= v.m_elem[0];
    m_elem[1] *= v.m_elem[1];
    m_elem[2] *= v.m_elem[2];
    m_elem[3] *= v.m_elem[3];

    return (*this);
}

template <typename m_type, uint16_t m_row>
inline CdtVector4<m_type, m_row>& CdtVector4<m_type, m_row>::operator /=(const CdtVector4& v)
{
    m_type den;

    den = v.m_elem[0];
    if (std::abs(den) < std::numeric_limits<m_type>::epsilon())
    {
        if (den < 0) den = -std::numeric_limits<m_type>::epsilon();
        else den = std::numeric_limits<m_type>::epsilon();
    }
    m_elem[0] /= den;

    den = v.m_elem[1];
    if (std::abs(den) < std::numeric_limits<m_type>::epsilon())
    {
        if (den < 0) den = -std::numeric_limits<m_type>::epsilon();
        else den = std::numeric_limits<m_type>::epsilon();
    }
    m_elem[1] /= den;

    den = v.m_elem[2];
    if (std::abs(den) < std::numeric_limits<m_type>::epsilon())
    {
        if (den < 0) den = -std::numeric_limits<m_type>::epsilon();
        else den = std::numeric_limits<m_type>::epsilon();
    }
    m_elem[2] /= den;

    den = v.m_elem[3];
    if (std::abs(den) < std::numeric_limits<m_type>::epsilon())
    {
        if (den < 0) den = -std::numeric_limits<m_type>::epsilon();
        else den = std::numeric_limits<m_type>::epsilon();
    }
    m_elem[3] /= den;

    return (*this);
}

template <typename m_type, uint16_t m_row>
inline CdtVector4<m_type, m_row>& CdtVector4<m_type, m_row>::operator =(const CdtVector<m_row, m_type>& v)
{
    m_elem[0] = v.m_elem[0];
    m_elem[1] = v.m_elem[1];
    m_elem[2] = v.m_elem[2];
    m_elem[3] = v.m_elem[3];

    return (*this);
}

template <typename m_type, uint16_t m_row>
inline CdtVector4<m_type, m_row>& CdtVector4<m_type, m_row>::operator +=(const CdtVector<m_row, m_type>& v)
{
    m_elem[0] += v.m_elem[0];
    m_elem[1] += v.m_elem[1];
    m_elem[2] += v.m_elem[2];
    m_elem[3] += v.m_elem[3];

    return (*this);
}

template <typename m_type, uint16_t m_row>
inline CdtVector4<m_type, m_row>& CdtVector4<m_type, m_row>::operator -=(const CdtVector<m_row, m_type>& v)
{
    m_elem[0] -= v.m_elem[0];
    m_elem[1] -= v.m_elem[1];
    m_elem[2] -= v.m_elem[2];
    m_elem[3] -= v.m_elem[3];

    return (*this);
}

template <typename m_type, uint16_t m_row>
inline CdtVector4<m_type, m_row>& CdtVector4<m_type, m_row>::operator *=(const CdtVector<m_row, m_type>& v)
{
    m_elem[0] *= v.m_elem[0];
    m_elem[1] *= v.m_elem[1];
    m_elem[2] *= v.m_elem[2];
    m_elem[3] *= v.m_elem[3];

    return (*this);
}

template <typename m_type, uint16_t m_row>
inline CdtVector4<m_type, m_row>& CdtVector4<m_type, m_row>::operator /=(const CdtVector<m_row, m_type>& v)
{
    m_type den;

    den = v.m_elem[0];
    if (std::abs(den) < std::numeric_limits<m_type>::epsilon())
    {
        if (den < 0) den = -std::numeric_limits<m_type>::epsilon();
        else den = std::numeric_limits<m_type>::epsilon();
    }
    m_elem[0] /= den;

    den = v.m_elem[1];
    if (std::abs(den) < std::numeric_limits<m_type>::epsilon())
    {
        if (den < 0) den = -std::numeric_limits<m_type>::epsilon();
        else den = std::numeric_limits<m_type>::epsilon();
    }
    m_elem[1] /= den;

    den = v.m_elem[2];
    if (std::abs(den) < std::numeric_limits<m_type>::epsilon())
    {
        if (den < 0) den = -std::numeric_limits<m_type>::epsilon();
        else den = std::numeric_limits<m_type>::epsilon();
    }
    m_elem[2] /= den;

    den = v.m_elem[3];
    if (std::abs(den) < std::numeric_limits<m_type>::epsilon())
    {
        if (den < 0) den = -std::numeric_limits<m_type>::epsilon();
        else den = std::numeric_limits<m_type>::epsilon();
    }
    m_elem[3] /= den;

    return (*this);
}

template<typename m_type, uint16_t m_row>
inline CdtVector4<m_type, m_row>& CdtVector4<m_type, m_row>::operator =(const CdtMatrix<m_row, 1, m_type>& v)
{
    m_elem[0] = v.m_elem[0];
    m_elem[1] = v.m_elem[1];
    m_elem[2] = v.m_elem[2];
    m_elem[3] = v.m_elem[3];

    return (*this);
}

template<typename m_type, uint16_t m_row>
inline CdtVector4<m_type, m_row>& CdtVector4<m_type, m_row>::operator +=(const CdtMatrix<m_row, 1, m_type>& v)
{
    m_elem[0] += v.m_elem[0];
    m_elem[1] += v.m_elem[1];
    m_elem[2] += v.m_elem[2];
    m_elem[3] += v.m_elem[3];

    return (*this);
}

template<typename m_type, uint16_t m_row>
inline CdtVector4<m_type, m_row>& CdtVector4<m_type, m_row>::operator -=(const CdtMatrix<m_row, 1, m_type>& v)
{
    m_elem[0] -= v.m_elem[0];
    m_elem[1] -= v.m_elem[1];
    m_elem[2] -= v.m_elem[2];
    m_elem[3] -= v.m_elem[3];

    return (*this);
}

template<typename m_type, uint16_t m_row>
inline CdtVector4<m_type, m_row>& CdtVector4<m_type, m_row>::operator *=(const CdtMatrix<m_row, 1, m_type>& v)
{
    m_elem[0] *= v.m_elem[0];
    m_elem[1] *= v.m_elem[1];
    m_elem[2] *= v.m_elem[2];
    m_elem[3] *= v.m_elem[3];

    return (*this);
}

template<typename m_type, uint16_t m_row>
inline CdtVector4<m_type, m_row>& CdtVector4<m_type, m_row>::operator /=(const CdtMatrix<m_row, 1, m_type>& v)
{
    m_type den;

    den = v.m_elem[0];
    if (std::abs(den) < std::numeric_limits<m_type>::epsilon())
    {
        if (den < 0) den = -std::numeric_limits<m_type>::epsilon();
        else den = std::numeric_limits<m_type>::epsilon();
    }
    m_elem[0] /= den;

    den = v.m_elem[1];
    if (std::abs(den) < std::numeric_limits<m_type>::epsilon())
    {
        if (den < 0) den = -std::numeric_limits<m_type>::epsilon();
        else den = std::numeric_limits<m_type>::epsilon();
    }
    m_elem[1] /= den;

    den = v.m_elem[2];
    if (std::abs(den) < std::numeric_limits<m_type>::epsilon())
    {
        if (den < 0) den = -std::numeric_limits<m_type>::epsilon();
        else den = std::numeric_limits<m_type>::epsilon();
    }
    m_elem[2] /= den;

    den = v.m_elem[3];
    if (std::abs(den) < std::numeric_limits<m_type>::epsilon())
    {
        if (den < 0) den = -std::numeric_limits<m_type>::epsilon();
        else den = std::numeric_limits<m_type>::epsilon();
    }
    m_elem[3] /= den;

    return (*this);
}

template<typename m_type, uint16_t m_row>
inline CdtVector4<m_type, m_row>& CdtVector4<m_type, m_row>::operator=(const m_type s)
{
    m_elem[0] = s;
    m_elem[1] = s;
    m_elem[2] = s;
    m_elem[3] = s;

    return (*this);
}

template <typename m_type, uint16_t m_row>
inline CdtVector4<m_type, m_row>& CdtVector4<m_type, m_row>::operator +=(const m_type s)
{
    m_elem[0] += s;
    m_elem[1] += s;
    m_elem[2] += s;
    m_elem[3] += s;

    return (*this);
}

template <typename m_type, uint16_t m_row>
inline CdtVector4<m_type, m_row>& CdtVector4<m_type, m_row>::operator -=(const m_type s)
{
    m_elem[0] -= s;
    m_elem[1] -= s;
    m_elem[2] -= s;
    m_elem[3] -= s;

    return (*this);
}

template <typename m_type, uint16_t m_row>
inline CdtVector4<m_type, m_row>& CdtVector4<m_type, m_row>::operator *=(const m_type s)
{
    m_elem[0] *= s;
    m_elem[1] *= s;
    m_elem[2] *= s;
    m_elem[3] *= s;

    return (*this);
}

template <typename m_type, uint16_t m_row>
inline CdtVector4<m_type, m_row>& CdtVector4<m_type, m_row>::operator /=(const m_type s)
{
    m_type den = s;

    if (std::abs(den) < std::numeric_limits<m_type>::epsilon())
    {
        if (den < 0) den = -std::numeric_limits<m_type>::epsilon();
        else den = std::numeric_limits<m_type>::epsilon();
    }

    m_elem[0] /= den;
    m_elem[1] /= den;
    m_elem[2] /= den;
    m_elem[3] /= den;

    return (*this);
}

template<typename m_type, uint16_t m_row>
inline CdtCommaInit<m_row, m_type> CdtVector4<m_type, m_row>::operator<<(const m_type s)
{
    m_elem[0] = s;
    return CdtCommaInit<m_row, m_type>(m_elem);
}

/* Arithmetic operators */
template <typename m_type, uint16_t m_row>
inline CdtVector4<m_type, m_row> CdtVector4<m_type, m_row>::operator -() const
{
    return CdtVector4(
        -m_elem[0],
        -m_elem[1],
        -m_elem[2],
        -m_elem[3]);
}

template <typename m_type, uint16_t m_row>
inline CdtVector4<m_type, m_row> CdtVector4<m_type, m_row>::operator +(const CdtVector4& v) const
{
    return CdtVector4(
        m_elem[0] + v.m_elem[0],
        m_elem[1] + v.m_elem[1],
        m_elem[2] + v.m_elem[2],
        m_elem[3] + v.m_elem[3]);
}

template <typename m_type, uint16_t m_row>
inline CdtVector4<m_type, m_row> CdtVector4<m_type, m_row>::operator -(const CdtVector4& v) const
{
    return CdtVector4(
        m_elem[0] - v.m_elem[0],
        m_elem[1] - v.m_elem[1],
        m_elem[2] - v.m_elem[2],
        m_elem[3] - v.m_elem[3]);
}

template <typename m_type, uint16_t m_row>
inline CdtVector4<m_type, m_row> CdtVector4<m_type, m_row>::operator *(const CdtVector4& v) const
{
    return CdtVector4(
        m_elem[0] * v.m_elem[0],
        m_elem[1] * v.m_elem[1],
        m_elem[2] * v.m_elem[2],
        m_elem[3] * v.m_elem[3]);
}

template <typename m_type, uint16_t m_row>
inline CdtVector4<m_type, m_row> CdtVector4<m_type, m_row>::operator /(const CdtVector4& v) const
{
    m_type den[4];

    den[0] = v.m_elem[0];
    den[1] = v.m_elem[1];
    den[2] = v.m_elem[2];
    den[3] = v.m_elem[3];

    if (std::abs(den[0]) < std::numeric_limits<m_type>::epsilon())
    {
        if (den[0] < 0) den[0] = -std::numeric_limits<m_type>::epsilon();
        else den[0] = std::numeric_limits<m_type>::epsilon();
    }
    if (std::abs(den[1]) < std::numeric_limits<m_type>::epsilon())
    {
        if (den[1] < 0) den[1] = -std::numeric_limits<m_type>::epsilon();
        else den[1] = std::numeric_limits<m_type>::epsilon();
    }
    if (std::abs(den[2]) < std::numeric_limits<m_type>::epsilon())
    {
        if (den[2] < 0) den[2] = -std::numeric_limits<m_type>::epsilon();
        else den[2] = std::numeric_limits<m_type>::epsilon();
    }
    if (std::abs(den[3]) < std::numeric_limits<m_type>::epsilon())
    {
        if (den[3] < 0) den[3] = -std::numeric_limits<m_type>::epsilon();
        else den[3] = std::numeric_limits<m_type>::epsilon();
    }

    return CdtVector4(
        m_elem[0] / den[0],
        m_elem[1] / den[1],
        m_elem[2] / den[2],
        m_elem[3] / den[3]);
}

template <typename m_type, uint16_t m_row>
inline CdtVector4<m_type, m_row> CdtVector4<m_type, m_row>::operator +(const CdtVector<m_row, m_type>& v) const
{
    return CdtVector4(
        m_elem[0] + v.m_elem[0],
        m_elem[1] + v.m_elem[1],
        m_elem[2] + v.m_elem[2],
        m_elem[3] + v.m_elem[3]);
}

template <typename m_type, uint16_t m_row>
inline CdtVector4<m_type, m_row> CdtVector4<m_type, m_row>::operator -(const CdtVector<m_row, m_type>& v) const
{
    return CdtVector4(
        m_elem[0] - v.m_elem[0],
        m_elem[1] - v.m_elem[1],
        m_elem[2] - v.m_elem[2],
        m_elem[3] - v.m_elem[3]);
}

template <typename m_type, uint16_t m_row>
inline CdtVector4<m_type, m_row> CdtVector4<m_type, m_row>::operator *(const CdtVector<m_row, m_type>& v) const
{
    return CdtVector4(
        m_elem[0] * v.m_elem[0],
        m_elem[1] * v.m_elem[1],
        m_elem[2] * v.m_elem[2],
        m_elem[3] * v.m_elem[3]);
}

template <typename m_type, uint16_t m_row>
inline CdtVector4<m_type, m_row> CdtVector4<m_type, m_row>::operator /(const CdtVector<m_row, m_type>& v) const
{
    m_type den[4];

    den[0] = v.m_elem[0];
    den[1] = v.m_elem[1];
    den[2] = v.m_elem[2];
    den[3] = v.m_elem[3];

    if (std::abs(den[0]) < std::numeric_limits<m_type>::epsilon())
    {
        if (den[0] < 0) den[0] = -std::numeric_limits<m_type>::epsilon();
        else den[0] = std::numeric_limits<m_type>::epsilon();
    }
    if (std::abs(den[1]) < std::numeric_limits<m_type>::epsilon())
    {
        if (den[1] < 0) den[1] = -std::numeric_limits<m_type>::epsilon();
        else den[1] = std::numeric_limits<m_type>::epsilon();
    }
    if (std::abs(den[2]) < std::numeric_limits<m_type>::epsilon())
    {
        if (den[2] < 0) den[2] = -std::numeric_limits<m_type>::epsilon();
        else den[2] = std::numeric_limits<m_type>::epsilon();
    }
    if (std::abs(den[3]) < std::numeric_limits<m_type>::epsilon())
    {
        if (den[3] < 0) den[3] = -std::numeric_limits<m_type>::epsilon();
        else den[3] = std::numeric_limits<m_type>::epsilon();
    }

    return CdtVector4(
        m_elem[0] / den[0],
        m_elem[1] / den[1],
        m_elem[2] / den[2],
        m_elem[3] / den[3]);
}

template<typename m_type, uint16_t m_row>
inline CdtVector4<m_type, m_row> CdtVector4<m_type, m_row>::operator +(const CdtMatrix<m_row, 1, m_type>& v) const
{
    return CdtVector4(
        m_elem[0] + v.m_elem[0],
        m_elem[1] + v.m_elem[1],
        m_elem[2] + v.m_elem[2],
        m_elem[3] + v.m_elem[3]);
}

template<typename m_type, uint16_t m_row>
inline CdtVector4<m_type, m_row> CdtVector4<m_type, m_row>::operator -(const CdtMatrix<m_row, 1, m_type>& v) const
{
    return CdtVector4(
        m_elem[0] - v.m_elem[0],
        m_elem[1] - v.m_elem[1],
        m_elem[2] - v.m_elem[2],
        m_elem[3] - v.m_elem[3]);
}

template<typename m_type, uint16_t m_row>
inline CdtVector4<m_type, m_row> CdtVector4<m_type, m_row>::operator *(const CdtMatrix<m_row, 1, m_type>& v) const
{
    return CdtVector4(
        m_elem[0] * v.m_elem[0],
        m_elem[1] * v.m_elem[1],
        m_elem[2] * v.m_elem[2],
        m_elem[3] * v.m_elem[3]);
}

template<typename m_type, uint16_t m_row>
inline CdtVector4<m_type, m_row> CdtVector4<m_type, m_row>::operator /(const CdtMatrix<m_row, 1, m_type>& v) const
{
    m_type den[4];

    den[0] = v.m_elem[0];
    den[1] = v.m_elem[1];
    den[2] = v.m_elem[2];
    den[3] = v.m_elem[3];

    if (std::abs(den[0]) < std::numeric_limits<m_type>::epsilon())
    {
        if (den[0] < 0) den[0] = -std::numeric_limits<m_type>::epsilon();
        else den[0] = std::numeric_limits<m_type>::epsilon();
    }
    if (std::abs(den[1]) < std::numeric_limits<m_type>::epsilon())
    {
        if (den[1] < 0) den[1] = -std::numeric_limits<m_type>::epsilon();
        else den[1] = std::numeric_limits<m_type>::epsilon();
    }
    if (std::abs(den[2]) < std::numeric_limits<m_type>::epsilon())
    {
        if (den[2] < 0) den[2] = -std::numeric_limits<m_type>::epsilon();
        else den[2] = std::numeric_limits<m_type>::epsilon();
    }
    if (std::abs(den[3]) < std::numeric_limits<m_type>::epsilon())
    {
        if (den[3] < 0) den[3] = -std::numeric_limits<m_type>::epsilon();
        else den[3] = std::numeric_limits<m_type>::epsilon();
    }

    return CdtVector4(
        m_elem[0] / den[0],
        m_elem[1] / den[1],
        m_elem[2] / den[2],
        m_elem[3] / den[3]);
}

template <typename m_type, uint16_t m_row>
inline CdtVector4<m_type, m_row> CdtVector4<m_type, m_row>::operator +(const m_type s) const
{
    return CdtVector4(
        m_elem[0] + s,
        m_elem[1] + s,
        m_elem[2] + s,
        m_elem[3] + s);
}

template <typename m_type, uint16_t m_row>
inline CdtVector4<m_type, m_row> CdtVector4<m_type, m_row>::operator -(const m_type s) const
{
    return CdtVector4(
        m_elem[0] - s,
        m_elem[1] - s,
        m_elem[2] - s,
        m_elem[3] - s);
}

template <typename m_type, uint16_t m_row>
inline CdtVector4<m_type, m_row> CdtVector4<m_type, m_row>::operator *(const m_type s) const
{
    return CdtVector4(
        m_elem[0] * s,
        m_elem[1] * s,
        m_elem[2] * s,
        m_elem[3] * s);
}

template <typename m_type, uint16_t m_row>
inline CdtVector4<m_type, m_row> CdtVector4<m_type, m_row>::operator /(const m_type s) const
{
    m_type den = s;

    if (std::abs(den) < std::numeric_limits<m_type>::epsilon())
    {
        if (den < 0) den = -std::numeric_limits<m_type>::epsilon();
        else den = std::numeric_limits<m_type>::epsilon();
    }

    return CdtVector4(
        m_elem[0] / den,
        m_elem[1] / den,
        m_elem[2] / den,
        m_elem[3] / den);
}

template <typename m_type, uint16_t m_row>
template <uint16_t col>
inline CdtMatrix<m_row, col, m_type> CdtVector4<m_type, m_row>::operator *(const CdtMatrix<1, col, m_type>& m) const
{
    m_type mat[m_row * col];
    uint16_t cnt;
    uint16_t irow, icol;

    for (irow = 0; irow < m_row; irow++)
    {
        for (cnt = col >> 2u, icol = 0; cnt > 0u; cnt--, icol += 4)
        {
            mat[irow * col + icol] = m_elem[irow] * m.m_elem[icol];
            mat[irow * col + icol + 1] = m_elem[irow] * m.m_elem[icol + 1];
            mat[irow * col + icol + 2] = m_elem[irow] * m.m_elem[icol + 2];
            mat[irow * col + icol + 3] = m_elem[irow] * m.m_elem[icol + 3];
        }

        for (cnt = col % 4u; cnt > 0u; cnt--, icol++)
            mat[irow * col + icol] = m_elem[irow] * m.m_elem[icol];
    }

    return CdtMatrix<m_row, col, m_type>(mat);
}

template <typename m_type, uint16_t m_row>
inline m_type CdtVector4<m_type, m_row>::dot(const CdtVector4& v) const
{
    return (
        m_elem[0] * v.m_elem[0] +
        m_elem[1] * v.m_elem[1] +
        m_elem[2] * v.m_elem[2] +
        m_elem[3] * v.m_elem[3]);
}

template <typename m_type, uint16_t m_row>
inline m_type CdtVector4<m_type, m_row>::dot(const CdtVector<m_row, m_type>& v) const
{
    return (
        m_elem[0] * v.m_elem[0] +
        m_elem[1] * v.m_elem[1] +
        m_elem[2] * v.m_elem[2] +
        m_elem[3] * v.m_elem[3]);
}

template<typename m_type, uint16_t m_row>
inline m_type CdtVector4<m_type, m_row>::dot(const CdtMatrix<m_row, 1, m_type>& v) const
{
    return (
        m_elem[0] * v.m_elem[0] +
        m_elem[1] * v.m_elem[1] +
        m_elem[2] * v.m_elem[2] +
        m_elem[3] * v.m_elem[3]);
}

/* Comparison operators */
template <typename m_type, uint16_t m_row>
inline bool CdtVector4<m_type, m_row>::operator ==(const CdtVector4& v) const
{
    if (std::abs(m_elem[0] - v.m_elem[0]) > m_tolerance) return false;
    if (std::abs(m_elem[1] - v.m_elem[1]) > m_tolerance) return false;
    if (std::abs(m_elem[2] - v.m_elem[2]) > m_tolerance) return false;
    if (std::abs(m_elem[3] - v.m_elem[3]) > m_tolerance) return false;

    return true;
}

template <typename m_type, uint16_t m_row>
inline bool CdtVector4<m_type, m_row>::operator !=(const CdtVector4& v) const
{
    if (std::abs(m_elem[0] - v.m_elem[0]) > m_tolerance) return true;
    if (std::abs(m_elem[1] - v.m_elem[1]) > m_tolerance) return true;
    if (std::abs(m_elem[2] - v.m_elem[2]) > m_tolerance) return true;
    if (std::abs(m_elem[3] - v.m_elem[3]) > m_tolerance) return true;

    return false;
}

template <typename m_type, uint16_t m_row>
inline bool CdtVector4<m_type, m_row>::operator ==(const CdtVector<m_row, m_type>& v) const
{
    if (std::abs(m_elem[0] - v.m_elem[0]) > m_tolerance) return false;
    if (std::abs(m_elem[1] - v.m_elem[1]) > m_tolerance) return false;
    if (std::abs(m_elem[2] - v.m_elem[2]) > m_tolerance) return false;
    if (std::abs(m_elem[3] - v.m_elem[3]) > m_tolerance) return false;

    return true;
}

template <typename m_type, uint16_t m_row>
inline bool CdtVector4<m_type, m_row>::operator !=(const CdtVector<m_row, m_type>& v) const
{
    if (std::abs(m_elem[0] - v.m_elem[0]) > m_tolerance) return true;
    if (std::abs(m_elem[1] - v.m_elem[1]) > m_tolerance) return true;
    if (std::abs(m_elem[2] - v.m_elem[2]) > m_tolerance) return true;
    if (std::abs(m_elem[3] - v.m_elem[3]) > m_tolerance) return true;

    return false;
}

template<typename m_type, uint16_t m_row>
inline bool CdtVector4<m_type, m_row>::operator==(const CdtMatrix<m_row, 1, m_type>& v) const
{
    if (std::abs(m_elem[0] - v.m_elem[0]) > m_tolerance) return false;
    if (std::abs(m_elem[1] - v.m_elem[1]) > m_tolerance) return false;
    if (std::abs(m_elem[2] - v.m_elem[2]) > m_tolerance) return false;
    if (std::abs(m_elem[3] - v.m_elem[3]) > m_tolerance) return false;

    return true;
}

template<typename m_type, uint16_t m_row>
inline bool CdtVector4<m_type, m_row>::operator!=(const CdtMatrix<m_row, 1, m_type>& v) const
{
    if (std::abs(m_elem[0] - v.m_elem[0]) > m_tolerance) return true;
    if (std::abs(m_elem[1] - v.m_elem[1]) > m_tolerance) return true;
    if (std::abs(m_elem[2] - v.m_elem[2]) > m_tolerance) return true;
    if (std::abs(m_elem[3] - v.m_elem[3]) > m_tolerance) return true;

    return false;
}

template <typename m_type, uint16_t m_row>
inline void CdtVector4<m_type, m_row>::Print(const char endChar)
{
#if defined(ARDUINO)
    for (uint16_t irow = 0; irow < m_row; irow++)
    {
        Serial.printf("%7.3f\n", (m_type)m_elem[irow]);
    }
    Serial.write(endChar);
#else
    for (uint16_t irow = 0; irow < m_row; irow++)
    {
        printf("%7.3f\n", (m_type)m_elem[irow]);
    }
    printf("%c", endChar);
#endif
}

//-- Template Function ------------------------------------------------------//
// scalar + vector
template <typename type, uint16_t row>
inline CdtVector4<type, row> operator +(const type s, const CdtVector4<type, row>& v)
{
    return CdtVector4<type, row>(
        v.m_elem[0] + s,
        v.m_elem[1] + s,
        v.m_elem[2] + s,
        v.m_elem[3] + s);
}

// scalar - vector
template <typename type, uint16_t row>
inline CdtVector4<type, row> operator -(const type s, const CdtVector4<type, row>& v)
{
    return CdtVector4<type, row>(
        s - v.m_elem[0],
        s - v.m_elem[1],
        s - v.m_elem[2],
        s - v.m_elem[3]);
}

// scalar * vector
template <typename type, uint16_t row>
inline CdtVector4<type, row> operator *(const type s, const CdtVector4<type, row>& v)
{
    return CdtVector4<type, row>(
        v.m_elem[0] * s,
        v.m_elem[1] * s,
        v.m_elem[2] * s,
        v.m_elem[3] * s);
}

// scalar / vector
template <typename type, uint16_t row>
inline CdtVector4<type, row> operator /(const type s, const CdtVector4<type, row>& v)
{
    type den[4];

    den[0] = v.m_elem[0];
    den[1] = v.m_elem[1];
    den[2] = v.m_elem[2];
    den[3] = v.m_elem[3];

    if (std::abs(den[0]) < std::numeric_limits<type>::epsilon())
    {
        if (den[0] < 0) den[0] = -std::numeric_limits<type>::epsilon();
        else den[0] = std::numeric_limits<type>::epsilon();
    }
    if (std::abs(den[1]) < std::numeric_limits<type>::epsilon())
    {
        if (den[1] < 0) den[1] = -std::numeric_limits<type>::epsilon();
        else den[1] = std::numeric_limits<type>::epsilon();
    }
    if (std::abs(den[2]) < std::numeric_limits<type>::epsilon())
    {
        if (den[2] < 0) den[2] = -std::numeric_limits<type>::epsilon();
        else den[2] = std::numeric_limits<type>::epsilon();
    }
    if (std::abs(den[3]) < std::numeric_limits<type>::epsilon())
    {
        if (den[3] < 0) den[3] = -std::numeric_limits<type>::epsilon();
        else den[3] = std::numeric_limits<type>::epsilon();
    }

    return CdtVector4<type, row>(
        s / den[0],
        s / den[1],
        s / den[2],
        s / den[3]);
}

typedef CdtVector4<> CdtVec4;
