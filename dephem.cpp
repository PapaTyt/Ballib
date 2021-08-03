//---------------------------------------------------------------------------
#pragma hdrstop
#include "dephem.h"

/* ����������� �� ���������. */
dph::EphemerisRelease::EphemerisRelease()
{
    /* ������������� ���������� ����������. */
    this->clear();
}

/* ����������� �� ���� � ��������� ����� ��������.
            ������ �����, �������� ���������� ��������. */
dph::EphemerisRelease::EphemerisRelease(const std::string& binaryFilePath)
{
    /* ������� ������� ���� ��������. */
    this->open(binaryFilePath);
}

/* ����������� �����������.
            �������� ���������� �������� � ������ � �����.
            ��� ��������� �������� ������ ���������. */
dph::EphemerisRelease::EphemerisRelease(const EphemerisRelease& other)
{
    if (other.m_ready)
    {
        copyHere(other);

        if (isDataCorrect())
        {
            m_ready = true;
        }
        else
        {
            m_ready = false;

            clear();
        }
    }
}

/* �������� �����������.
            �������� ���������� �������� � ������ � �����.
            ��� ��������� �������� ������ ���������. */
dph::EphemerisRelease&
dph::EphemerisRelease::operator=(const EphemerisRelease& other)
{
    if (other.m_ready)
    {
        clear();
        copyHere(other);

        if (isDataCorrect())
        {
            m_ready = true;
        }
        else
        {
            m_ready = false;

            clear();
        }
    }

    return *this;
}

/* ����������.
            ������ ����������.*/
dph::EphemerisRelease::~EphemerisRelease()
{
    m_binaryFileStream.close();
}

/* ������� ������� ���� ��������. */
void dph::EphemerisRelease::open(const std::string& binaryFilePath)
{
    /* ������������� ���������� ����������. */
    clear();

    /* ����������� ���� � �����. */
    m_binaryFilePath = binaryFilePath;

    /* �������� �����. */
    m_binaryFileStream.open(m_binaryFilePath.c_str(), std::ios::binary);

    if (m_binaryFileStream.is_open())
    {
        readAndPackData();

        if (isDataCorrect())
        {
            m_ready = true;
        }
        else
        {
            clear();
        }
    }
    else
    {
        m_binaryFilePath.clear();
    }
}

/* ������-������ (��� ������ ���������) ������ ���� ������������ �������. */
void dph::EphemerisRelease::calculateBody(unsigned calculationResult,
	unsigned targetBody, unsigned centerBody, double JED,
	double* resultArray) const
{
    /*
     *  ���������:
     *   calculationResult � ������ ���������� ����������.
     *   targetBody        � ���������� ����� �������� ����.
     *   centerBody        � ���������� ����� ������������ ����.
     *   JED               � ������ ������� (Julian Epehemris Date).
     *   resultArray       � ��������� ����������.
     *
     *  ����������:
     *    � ���� � ����� ������ �������� ���������, �� �� ������ ��������.
     */

    /*
     *  ���������� �������� ����������
     *
     *  calculationResult:
     *		1 - �������� �������� ������-�������,
     *		2 - �������� �������� ������� ���������
     *
     *   targetBody, centerBody:
     *      1  � ��������;
     *	    2  � ������;
     *		3  � �����;
     *		4  � ����;
     *		5  � ������;
     *		6  � ������;
     *		7  � ����;
     *		8  � ������;
     *		9  � ������;
     *		10 � ����;
     *		11 � ������;
     *		12 � ��������� ��������� �������;
     *		13 � ��������� ������� �����-����;
     *
     *  JED:
     *      JED ������ ������������ ����������: [m_startDate : m_endDate].
     *
     *  resultArray:
     *      � ������ �������:
     *          3 ��� calculationResult = 0;
     *          6 ��� calculationResult = 1.
     *      � �� ������ ���� ������� ����������.
     */

    /* �������� ���������� ����������. */
    if (this->m_ready == false)
    {
        return;
    }
    else if (calculationResult > 1)
    {
        return;
    }
    else if (targetBody == 0 || centerBody == 0)
    {
        return;
    }
    else if (targetBody > 13 || centerBody > 13)
    {
        return;
    }
    else if (JED < m_startDate || JED > m_endDate)
    {
        return;
    }
	else if (resultArray == NULL)
    {
        return;
	}

    /* ���������� ��������� ���������. */
    unsigned componentsCount = (calculationResult == CALC_STATE) ? 6 : 3;

    /* ����� �������� ���������� � ����������� �� ���������� �������� �
       ������������ ����. */
    if (targetBody == centerBody)
    {
        /* ������ #1: ������� ���� �������� �����������.
                      ����������� �������� ������� ������.*/

        /* ���������� ������� ������ */
        std::memset(resultArray, 0, sizeof(double) * componentsCount);
    }
    else if (targetBody == B_SSBARY || centerBody == B_SSBARY)
    {
        /* ������ #2: ������� ��� ����������� ����� �������� ��������� ��.
                      ��� ��� ��� ������ calculateBase ��� ��� ���������� ������
                      ������������ ��������� ��, �� ���������� ���������
                      ���������� ������� ����. � ������, ���� ������� �����
                      �������� ��� ��������� ��, �� ������������ "����������"
                      ������ ������� ����. */

        /* ������ ����, ��� �� �������� ����������� ��. */
        unsigned notSSBARY = targetBody == B_SSBARY ? centerBody : targetBody;

        /* ����� ������ ���������� � ����������� �� ����. */
        switch (notSSBARY)
        {
            case B_EARTH:
                calculateBaseEarth(JED, calculationResult, resultArray);
                break;

            case B_MOON:
                calculateBaseMoon(JED, calculationResult, resultArray);
                break;

            case B_EMBARY:
                calculateBaseItem(2, JED, calculationResult, resultArray);
                break;

            default:
                calculateBaseItem(notSSBARY - 1, JED, calculationResult,
                    resultArray);
        }

        /* ���� ��������� �� �������� ������� �����, �� ������������
           "����������" ������. */
        if (targetBody == B_SSBARY)
        {
            for (unsigned i = 0; i < componentsCount; ++i)
            {
                resultArray[i] = -resultArray[i];
            }
        }
    }
    else if (targetBody * centerBody == 30 && targetBody + centerBody == 13)
    {
        /* ������ #3: ������� � ����������� ������ �������� ����� � ����.
                      � ���� ������ ���������� �������� �������� ��������� ����
                      ������������ ����� (������� ������� #9 (�� ����).
                      � ������, ���� ������� ����� �������� �����, ��
                      ������������ "���������� ������" */

        /* ��������� ������-������� (��� ������� ���������) ���� ������������
           �����. */
        calculateBaseItem(9, JED, calculationResult, resultArray);

        /* ���� ������� ����� �������� �����, �� ������������ "����������"
           ������. */
        if (targetBody == B_EARTH)
        {
            for (unsigned i = 0; i < componentsCount; ++i)
            {
                resultArray[i] = -resultArray[i];
            }
        }
    }
    else
    {
        /* ������ #4: ��� ��������� ���������� ���.
                      ������� ����������� �������� ������������ ����
                      ������������ ���������� ��, ����� - ��������. �����������
                      �������� ������� ����� �������� ������������ ���� �
                      ��������. */

        /* ������ ��� ������������ ����. */
        double centerBodyArray[6];

        /* ��� ��������. */
        for (unsigned i = 0; i <= 1; ++i)
        {
            /* ����������� ������� � ������� � ����������� �� ������ ��������.
               i == 0 : ������ � ����������� �����.
               i == 1 : ������ � ������� �����. */
            unsigned currentBodyIndex = i == 0 ? centerBody : targetBody;
            double* currentArray = i == 0 ? centerBodyArray : resultArray;

            /* ����� ������ ���������� � ����������� �� ����. */
            switch (currentBodyIndex)
            {
                case B_EARTH:
                    calculateBaseEarth(JED, calculationResult, currentArray);
                    break;

                case B_MOON:
                    calculateBaseMoon(JED, calculationResult, currentArray);
                    break;

                case B_EMBARY:
                    calculateBaseItem(2, JED, calculationResult, currentArray);
                    break;

                default:
                    calculateBaseItem(currentBodyIndex - 1, JED,
                        calculationResult, currentArray);
            }
        }

        /* ������� ����� �������� ������������ � �������� ����. */
        for (unsigned i = 0; i < componentsCount; ++i)
        {
            resultArray[i] -= centerBodyArray[i];
        }
    }
}
/* ������-������ (��� ������ ���������) ������ ���� ������������ �������. */
void dph::EphemerisRelease::calculateR(unsigned targetBody, unsigned centerBody,
					double JED, double r[3]) const
{

calculateBody(1, targetBody, centerBody, JED, r);

}
/* ������-������ (��� ������ ���������) ������ ���� ������������ �������. */
void dph::EphemerisRelease::calculateRV(unsigned targetBody, unsigned centerBody,
					double JED, double r[3], double v[3]) const
{
double rv[6];
calculateBody(2, targetBody, centerBody, JED, rv);
for(int i=0; i<3; i++){
	r[i]=rv[i];
	v[i]=rv[i+3];
}
}



/* �������� ������ �� �������������� ���������, ���������� � ������� DE. */
void dph::EphemerisRelease::calculateOther(unsigned calculationResult,
    unsigned otherItem, double JED, double* resultArray) const
{
    /*
     *  ���������:
     *   calculationResult � ������ ���������� ����������.
     *   otherItem         � ���������� ����� �������� ��������
     *   JED               � ������ ������� (Julian Epehemris Date).
     *   resultArray       � ��������� ����������.
     *
     *  ����������:
     *    � ���� � ����� ������ �������� ���������, �� �� ������ ��������.
     *    � �� ������ � ������� �������� �������� ������������� ����, ���������
     *      � ��� ������� ����� ��������.
     */

    /*
     *  ���������� �������� ����������
     *
     *   calculationResult:
     *		1 - �������� �������� ������-�������,
     *		2 - �������� �������� ������� ���������
     *
     *   otherItem:
     *      14 � ������ ������� � �������� � ������� (������ IAU 1980);
     *      15 � �������� ������ ������;
     *		16 � ������� �������� ������ ������;
     *		17 � TT-TDB (� ������ �����).
     *
     *  JED:
     *      JED ������ ������������ ����������: [m_startDate : m_endDate].
     *
     *  resultArray:
     *      � ������ ������� � ����������� �� ���������� ���������� � ��������
     *        ��������.
     *      � �� ������ ���� ������� ����������.
     */

    /* �������� ���������� ����������. */
    if (this->m_ready == false)
    {
        return;
    }
    else if (calculationResult > 1)
    {
        return;
    }
    else if (otherItem < 14 || otherItem > 17)
    {
        return;
    }
    else if (JED < m_startDate || JED > m_endDate)
    {
        return;
    }
	else if (resultArray == NULL)
    {
        return;
    }
    else
    {
        calculateBaseItem(otherItem - 3, JED, calculationResult, resultArray);
    }
}


/* ���������� ������� � �������������. */
bool dph::EphemerisRelease::isReady() const
{
    return m_ready;
}

/* ������ ��������� ���� ��� ���������. */
double dph::EphemerisRelease::startDate() const
{
    return m_startDate;
}

/* ��������� ��������� ���� ��� ���������. */
double dph::EphemerisRelease::endDate() const
{
    return m_endDate;
}

/* ����� �������. */
uint32_t dph::EphemerisRelease::releaseIndex() const
{
    return m_releaseIndex;
}

/* ������������ ���������� DE. */
const std::string& dph::EphemerisRelease::releaseLabel() const
{
    return m_releaseLabel;
}

/* �������� ��������� �� � �����. */
double dph::EphemerisRelease::constant(const std::string& constantName) const
{
    if (m_ready == false)
    {
        return 0.0;
    }
    else if (constantName == "AU")
    {
        return m_au;
    }
    else if (constantName == "EMRAT")
    {
        return m_emrat;
    }
    else if (constantName == "DENUM")
    {
        return m_releaseIndex;
    }
    else
    {
		std::map<std::string, double>::const_iterator found =
			m_constants.find(constantName);

		if (found != m_constants.end())
		{
			return found->second;
		}
		else
		{
			return 0;
        }
    }
}

/* �������� ������������� ������� (' ') � ����� ������� �������� "charArray"
       ������� "arraySize". */
std::string dph::EphemerisRelease::cutBackSpaces(const char* charArray,
    size_t arraySize)
{
    for (size_t i = arraySize - 1; i > 0; --i)
    {
        if (charArray[i] == ' ' && charArray[i - 1] != ' ')
        {
            return std::string(charArray, i);
        }
    }

    return std::string(charArray, arraySize);
}

/* ���������� ������� � ������������ ���������. */
void dph::EphemerisRelease::clear()
{
    m_ready = false;

    m_binaryFilePath.clear();
    m_binaryFileStream.close();

    m_releaseLabel.clear();
    m_releaseIndex = 0;
    m_startDate = 0.0;
    m_endDate = 0.0;
    m_blockTimeSpan = 0.0;
    std::memset(m_keys, 0, sizeof(m_keys));
    m_au = 0.0;
    m_emrat = 0.0;
    std::map<std::string, double>().swap(m_constants);

    m_blocksCount = 0;
    m_ncoeff = 0;
    m_dimensionFit = 0;
    m_blockSize_bytes = 0;

    std::vector<double>().swap(m_buffer);
    std::vector<double>(1).swap(m_poly);
    std::vector<double>(2).swap(m_dpoly);

    m_poly[0]  = 1;
    m_dpoly[0] = 0;
    m_dpoly[1] = 1;
}

/* ����������� ���������� �� ������� "other" � ������� ������. */
void dph::EphemerisRelease::copyHere(const EphemerisRelease& other)
{
    /* ������������ �:
        - ������������ �����������.
        - ��������� �����������. */

    m_ready = other.m_ready;

    m_binaryFilePath	= other.m_binaryFilePath;

    m_binaryFileStream.close();
    m_binaryFileStream.open(other.m_binaryFilePath.c_str(), std::ios::binary);

    m_releaseLabel =	other.m_releaseLabel;
    m_releaseIndex =	other.m_releaseIndex;
    m_startDate =		other.m_startDate;
    m_endDate =			other.m_endDate;
    m_blockTimeSpan =	other.m_blockTimeSpan;
    std::memcpy(m_keys, other.m_keys, sizeof(m_keys));
    m_au =				other.m_au;
    m_emrat =			other.m_emrat;
    m_constants =		other.m_constants;

    m_blocksCount =		other.m_blocksCount;
    m_ncoeff =			other.m_ncoeff;
    m_emrat2 =			other.m_emrat2;
    m_dimensionFit =	other.m_dimensionFit;
    m_blockSize_bytes = other.m_blockSize_bytes;

    m_buffer =	other.m_buffer;
    m_poly =	other.m_poly;
    m_dpoly =	other.m_poly;
}

/* ������� ������ �� ��������� ����� ������� DE � ������. */
void dph::EphemerisRelease::readAndPackData()
{
    /* ��������� ���������� � ������� DE (�����). */
	char releaseLabel_buffer[DPH_RLS_LABELS_COUNT][DPH_RLS_LABEL_SIZE];

	/* ����� ��������. */
	char constantsNames_buffer[DPH_CCOUNT_MAX_NEW][DPH_CNAME_SIZE];

	/* �������� ��������. */
	double constantsValues_buffer[DPH_CCOUNT_MAX_NEW];

    /* ���������� �������� � ����� ��������. */
	uint32_t constantsCount;


    /* [������ �����] */

    m_binaryFileStream.seekg(0, std::ios::beg);
    m_binaryFileStream.read(reinterpret_cast<char*>(&releaseLabel_buffer),
							DPH_RLS_LABEL_SIZE * DPH_RLS_LABELS_COUNT);
    m_binaryFileStream.read(reinterpret_cast<char*>(&constantsNames_buffer),
							DPH_CNAME_SIZE * DPH_CCOUNT_MAX_OLD);
    m_binaryFileStream.read(reinterpret_cast<char*>(&m_startDate), 8);
    m_binaryFileStream.read(reinterpret_cast<char*>(&m_endDate), 8);
    m_binaryFileStream.read(reinterpret_cast<char*>(&m_blockTimeSpan), 8);
    m_binaryFileStream.read(reinterpret_cast<char*>(&constantsCount), 4);
    m_binaryFileStream.read(reinterpret_cast<char*>(&m_au), 8);
    m_binaryFileStream.read(reinterpret_cast<char*>(&m_emrat), 8);
    m_binaryFileStream.read(reinterpret_cast<char*>(&m_keys), (12 * 3) * 4);
    m_binaryFileStream.read(reinterpret_cast<char*>(&m_releaseIndex), 4);
    m_binaryFileStream.read(reinterpret_cast<char*>(&m_keys[12]), (3) * 4);

    /* ������ �������������� ��������. */
    if (constantsCount > 400)
    {
        /* ���������� �������������� ��������. */
		int extraConstantsCount = constantsCount - DPH_CCOUNT_MAX_OLD;
		int bufsize = extraConstantsCount * static_cast<int>(DPH_CNAME_SIZE);

        m_binaryFileStream.read(
			reinterpret_cast<char*>(&constantsNames_buffer[DPH_CCOUNT_MAX_OLD]),
            bufsize);
    }

    /* ������ �������������� ������. */
    m_binaryFileStream.read(reinterpret_cast<char*>(&m_keys[13]), (3 * 2) * 4);

    /* ������� ncoeff (���������� ������������� � �����). */
    m_ncoeff = 2;
    for (int i = 0; i < 15; ++i)
    {
        /* ���������� ��������� ��� ���������� ��������. */
        int comp = i == 11 ? 2 : i == 14 ? 1 : 3;
        m_ncoeff += static_cast<uint32_t>(comp) * m_keys[i][1] * m_keys[i][2];
    }

    /* ������� � ����� � ����������� � �� ������. */
	if (constantsCount <= DPH_CCOUNT_MAX_NEW)
    {
        int offset  = static_cast<int>(m_ncoeff * 8);
        int bufsize = static_cast<int>(constantsCount * 8);
        m_binaryFileStream.seekg(offset, std::ios::beg);
        m_binaryFileStream.read(
            reinterpret_cast<char*>(&constantsValues_buffer), bufsize);
    }


    /* [�������������� � �������� ��������� ����������] */

    /* ������������ ����� ����� ������� DE. */
	for (size_t i = 0; i < DPH_RLS_LABELS_COUNT; ++i)
    {
		m_releaseLabel += cutBackSpaces(releaseLabel_buffer[i],
			DPH_RLS_LABEL_SIZE);
        m_releaseLabel += '\n';
    }

    /* ���������� ���������� m_constants ������� � ���������� ��������. */
	if (constantsCount > 0 && constantsCount <= DPH_CCOUNT_MAX_NEW)
    {
        for (uint32_t i = 0; i < constantsCount; ++i)
        {
            std::string constantName = cutBackSpaces(constantsNames_buffer[i],
				DPH_CNAME_SIZE);
            m_constants[constantName] = constantsValues_buffer[i];
        }
    }

    /* �������������� ����������. */
    additionalCalculations();
}

/* �������������� ���������� ����� ������ �����.
   ������ � ������ readAndPackData(). */
void dph::EphemerisRelease::additionalCalculations()
{
    /* �������������� ������������. */
    m_emrat2 = 1.0 / (1 + m_emrat);
    m_dimensionFit = 1.0 / (43200.0 * m_blockTimeSpan);

    /* ���������� ������ � ��������������. */
    m_blocksCount = size_t((m_endDate - m_startDate) / m_blockTimeSpan);

    /* ������������ ������� ����������������� ��������. */
    size_t maxPolynomsCount = 0;
    for (int i = 0; i < 15; ++i)
    {
        if (m_keys[i][1] > maxPolynomsCount)
        {
            maxPolynomsCount = m_keys[i][1];
        }
    }

    /* ������ ����� � ������. */
    m_blockSize_bytes = m_ncoeff * sizeof(double);

    /* �������������� ������. */
    m_buffer.resize(m_ncoeff);
    m_poly.resize(maxPolynomsCount);
    m_dpoly.resize(maxPolynomsCount);
}

/* �������� ��������, ���������� � ������� � �������� �����. */
bool dph::EphemerisRelease::isDataCorrect() const
{
    /* � ������ ������ ����������� ������ �� ���������, ������� ����� ��������
       ��������������� �� ���������� �������� ���������, ���������� � �������
       ��������. */

    if (m_binaryFileStream.is_open() == false)       return false;
    if (m_startDate >= m_endDate)                    return false;
    if (m_blockTimeSpan == 0)                        return false;
    if ((m_endDate - m_startDate) < m_blockTimeSpan) return false;
    if (m_emrat == 0)                                return false;
    if (m_ncoeff == 0)                               return false;
    if (check_blocksDates() == false)                return false;

    return true;
}


/* �������� ��������� � �������� ��� ���� ������ � �����.
   ������������ ����������� ����� � ����������� ���� �������������.
   ������ � ������ �������� isDataCorrect(). */
bool dph::EphemerisRelease::check_blocksDates() const
{
    /* ����� ������� ����� � �������������� � �����. */
    size_t firstBlockAdress = m_blockSize_bytes * 2;

    /* ������� � ������� �����. */
    m_binaryFileStream.seekg(static_cast<int>(firstBlockAdress), std::ios::beg);

    /* �������� ����� ������� ����� ������ ���� ������ �������������. */
    size_t subBlockOffset = (m_ncoeff - 2) * sizeof(double);

    for (size_t blockIndex = 0; blockIndex < m_blocksCount; ++blockIndex)
    {
        /* ������ ��� ������ ������ ���� ������������� �� �������� �����. */
        double blockDates[2] = {0.0, 0.0};

        /* ������. */
        m_binaryFileStream.read(reinterpret_cast<char*>(& blockDates),
            sizeof(blockDates));

        /* ��������, ������� ������ ����. */
        double blockStartDate = m_startDate + blockIndex * m_blockTimeSpan;
        double blockEndDate = blockStartDate + m_blockTimeSpan;

        if (blockDates[0] != blockStartDate || blockDates[1] != blockEndDate)
        {
            return false;
        }

        /* ������� � ���������� �����. */
        m_binaryFileStream.seekg(static_cast<int>(subBlockOffset),
            std::ios::cur);
    }

    return true;
}

/* ���������� ������� "m_buffer" �������������� ���������� �����. */
void dph::EphemerisRelease::fillBuffer(size_t block_num) const
{
    size_t adress  = (2 + block_num) * m_blockSize_bytes;
    int    bufsize = static_cast<int>(m_ncoeff * 8);

    m_binaryFileStream.seekg(static_cast<int>(adress), std::ios::beg);
    m_binaryFileStream.read(reinterpret_cast<char*>(&m_buffer[0]),
        bufsize);
}

/* ������������ ��������� �������� ��������. */
void dph::EphemerisRelease::interpolatePosition(unsigned baseItemIndex,
    double normalizedTime, const double* coeffArray, unsigned componentsCount,
    double* resultArray) const
{
    /* ���������� ������������� �� ����������. */
    uint32_t cpec = m_keys[baseItemIndex][1];

    /* ��������������� ���������� ��������� (���������� �� ����). */
    m_poly[1] = normalizedTime;

    /* ���������� ��������� (���������� �� ����). */
    for (uint32_t i = 2; i < cpec; ++i)
    {
        m_poly[i] = 2 * normalizedTime * m_poly[i - 1] - m_poly[i - 2];
    }

    /* ��������� ������� ���������� ����������. */
    memset(resultArray, 0, sizeof(double) * componentsCount);

    /* ���������� ���������. */
    for (unsigned i = 0; i < componentsCount; ++i)
    {
        for (uint32_t j = 0; j < cpec; ++j)
        {
            resultArray[i] += m_poly[j] * coeffArray[i * cpec + j];
        }
    }
}

/* ������������ ��������� � �� ������ ����������� �������� ��������. */
void dph::EphemerisRelease::interpolateState(unsigned baseItemIndex,
    double normalizedTime, const double* coeffArray, unsigned componentsCount,
    double* resultArray) const
{
    /* ���������� ������������� �� ����������. */
    uint32_t cpec = m_keys[baseItemIndex][1];

    /* ��������������� ���������� ��������� (���������� �� ����). */
    m_poly[1]  = normalizedTime;
    m_poly[2]  = 2 * normalizedTime * normalizedTime - 1;
    m_dpoly[2] = 4 * normalizedTime;

    /* ���������� ��������� (���������� �� ����). */
    for (uint32_t i = 3; i < cpec; ++i)
    {
        m_poly[i] = 2 * normalizedTime *  m_poly[i - 1] -  m_poly[i - 2];
        m_dpoly[i] = 2 * m_poly[i - 1] + 2 * normalizedTime * m_dpoly[i - 1] -
            m_dpoly[i - 2];
    }

    /* ��������� ������� ���������� ����������. */
    memset(resultArray, 0, sizeof(double) * componentsCount * 2);

    /* ����������� ��� ���������� �����������. */
    double derivative_units = m_keys[baseItemIndex][2] * m_dimensionFit;

    /* ���������� ���������. */
    for (unsigned i = 0; i < componentsCount; ++i)
    {
        for (uint32_t j = 0; j < cpec; ++j, ++coeffArray)
        {
            resultArray[i]                   +=  m_poly[j] * *coeffArray;
            resultArray[i + componentsCount] += m_dpoly[j] * *coeffArray;
        }

        resultArray[i + componentsCount] *= derivative_units;
    }
}

/* �������� ��������� �������� ��������. */
void dph::EphemerisRelease::calculateBaseItem(unsigned baseItemIndex,
    double JED, unsigned calculationResult, double* resultArray) const
{
    /*
     *  ���������:
     *    baseItemIndex     � ������ �������� �������� ������� DE (�� ����).
     *    JED               � ������ ������� (Julian Epehemris Date).
     *    calculationResult � ������ ���������� ����������.
     *    resultArray       � ��������� ����������.
     */

    /*
     *  ���������� �������� ����������
     *
     *  baseItemIndex:
     *      0  � ��������;
     *      1  � ������;
     *      2  � ��������� ������� �����-����;
     *      3  � ����;
     *      4  � ������;
     *      5  � ������;
     *      6  � ����;
     *      7  � ������;
     *      8  � ������;
     *      9  � ���� (������������ �����);
     *      10 � ������;
     *      11 � ������ ������� � �������� � ������� (������ IAU 1980);
     *      12 � �������� ������ ������;
     *      13 � ������� �������� ������ ������;
     *      14 � TT-TDB (� ������ �����).
     *
     *  JED:
     *      JED ������ ������������ ����������: [m_startDate : m_endDate].
     *
     *  resultArray:
     *      � ������ ������� � ����������� �� ���������� ���������� � ��������
     *        ��������.
     *      � �� ������ ���� ������� ����������.
     */

    /* ��������!
       � ���� ���������� ������� ����� ���������� "normalizedTime" � "offset"
       ����� ��������. */

    /* ����. ����� ������������ ���� ������ � �������. */
    double normalizedTime = (JED - m_startDate) / m_blockTimeSpan;

    /* ���������� ����� �����, �����. �������� ���� JED
       (����� ����� �� normalizedTime). */
    size_t offset = static_cast<size_t>(normalizedTime);

    /* ���������� ������� �������������� ���������� �����.
        ���� ��������� ���� ��� � ���� �������, �� �� �� ����������� ��������.
        m_buffer[0] - ���� ������ �����.
        m_buffer[1] - ���� ��������� �����. */
    if (JED < m_buffer[0] || JED >= m_buffer[1])
    {
        /* ���� JED ����� ��������� ���������� ���� ��� ����������, ��
           ����������� ��������� ����. */

        fillBuffer(offset - (JED == m_endDate ? 1 : 0));
    }

    if (JED == m_endDate)
    {
        /* ���������� ����� �������� (��������� �������). */
        offset = m_keys[baseItemIndex][2] - 1;

        /* ����. ����� ������������ �������� (� ��������� �� -1 �� 1). */
        normalizedTime = 1;
    }
    else
    {
        /* ����. ����� ������������ ���� ���������. */
        normalizedTime = (normalizedTime - offset) * m_keys[baseItemIndex][2];

        /* ���������� ����� �������� (����� ����� �� normalizedTime). */
        offset = static_cast<size_t>(normalizedTime);

        /* ����. ����� ������������ �������� (� ��������� �� -1 �� 1). */
        normalizedTime = 2 * (normalizedTime - offset) - 1;
    }

    /* ���������� ��������� ��� ���������� �������� ��������. */
    unsigned componentsCount = (baseItemIndex == 11) ? 2 :
                               (baseItemIndex == 14) ? 1 : 3;

    /* ���������� ����� ������� ������������ � �����. */
    size_t coeff_pos = m_keys[baseItemIndex][0] - 1 +
        componentsCount * offset * m_keys[baseItemIndex][1];

    /* ����� ������ ���������� � ����������� �� ��������� ����������
       ����������. */
    switch(calculationResult)
    {
    case CALC_POS :
        interpolatePosition(baseItemIndex, normalizedTime, &m_buffer[coeff_pos],
            componentsCount, resultArray);
        break;

    case CALC_STATE :
        interpolateState(baseItemIndex, normalizedTime, &m_buffer[coeff_pos],
            componentsCount, resultArray);
        break;

    default:
        memset(resultArray, 0, componentsCount * sizeof(double));
    }
}

/* �������� �������� ������-������� (��� ������� ���������) �����
   ������������ ���������� ��������� �������. */
void dph::EphemerisRelease::calculateBaseEarth(double JED, unsigned
    calculationResult, double* resultArray) const
{
    /* ������-������ (��� ������ ���������) ���������� �������� �����-����
       ������������ ���������� ��. */
    calculateBaseItem(2, JED, calculationResult, resultArray);

    /* ������-������ (��� ������ ���������) ���� ����������� �����. */
    double MoonRelativeEarth[6];
    calculateBaseItem(9, JED, calculationResult, MoonRelativeEarth);

    /* ���������� ���������. */
    unsigned componentsCount = calculationResult == CALC_POS ? 3 : 6;

    /* ������-������ (��� ������ ���������) ����� ������������ ���������� ��. */
    for (unsigned i = 0; i < componentsCount; ++i)
    {
        resultArray[i] -= MoonRelativeEarth[i] * m_emrat2;
    }
}

/* ������-������ (��� ������ ���������) ���� ������������ ���������� ��. */
void dph::EphemerisRelease::calculateBaseMoon(double JED, unsigned
    calculationResult, double* resultArray) const
{
    /* ������-������ (��� ������ ���������) ���������� �������� �����-����
       ������������ ���������� ��. */
    calculateBaseItem(2, JED, calculationResult, resultArray);

    /* ������-������ (��� ������ ���������) ���� ����������� �����. */
    double MoonRelativeEarth[6];
    calculateBaseItem(9, JED, calculationResult, MoonRelativeEarth);

    /* ���������� ���������. */
    unsigned componentsCount = calculationResult == CALC_POS ? 3 : 6;

    /* ������������� ���������. */
    for (unsigned i = 0; i < componentsCount; ++i)
    {
        resultArray[i] += MoonRelativeEarth[i] * (1 - m_emrat2);
    }
}

//---------------------------------------------------------------------------

#pragma package(smart_init)
