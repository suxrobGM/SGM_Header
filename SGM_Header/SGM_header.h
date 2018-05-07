/**
 *  Useful functions
 *  In this header used some standart functions of C/C++ also used C++11 features
 *  This header created by SuxrobGM
 *
 *  Functions:
 *  1) Factorial()
 *  2) CountDigits()
 *  3) DivideNumsByDig()    (+2 overloaded)
 *  4) GetMax()             (+3 overloaded)
 *  5) GetMin()             (+3 overloaded)
 *  6) SortInc_B()
 *  7) SortDec_B()
 *  8) PrintArray()
 *  9) ScanArray()
 *  10) Recursive_BinarySearch()
 *  11) BinarySearch()
 *  12) GCD()
 *  13) LCM()
 *
 *  Classes:
 *  1) HashEntry
 *	2) HashTable
 *  3) Point
 *  4) Determinant      (Not completed)
 *  5) BigInteger       (Not completed)
 *  6) Stack
 *  7) Queue
 *  8) Deque
 *  9) Node
 *  10) ForwardList
 *  11) BinaryTree
 *  12) Matrix
 *  13) Graph ?
 *
 */
#ifndef SGM_HEADER_H_INCLUDED
#define SGM_HEADER_H_INCLUDED
#define DEBUG 1

#include <iostream>
#include <iomanip>
#include <string>
#include <math.h>
#include <stdlib.h>

namespace SGM {

    typedef unsigned long long int UINT64;
    typedef unsigned long int UINT32;
    typedef unsigned int UINT16;
    typedef unsigned short int USHORT;
    typedef long long int INT64;
    typedef long int INT32;
    typedef int INT16;
    typedef short int SHORT;
    typedef float DWORD;
    typedef double REAL;
    typedef long double LREAL;
    typedef unsigned char BYTE;
    typedef char CHAR;
    typedef wchar_t WCHAR;
    typedef size_t UINT_T;

//-------------------------------------------------------------------------------------------------------------------------------
//----------------------------------------------  FUNCTIONS  --------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------------------------
    /** Factorial  n! */
    UINT64 Factorial(UINT_T InputNum)
    {
        UINT64 result=1;
        for(UINT32 i=1; i<=InputNum; i++)
        {
            result*=i;
        }
        return result;
    }

    /** Count digits in the integer number */
    template<typename INT_T>
    UINT64 CountDigits(INT_T inputNum)
    {
        UINT32 Count=0;
        while(inputNum !=0)
        {
            inputNum /= 10;
            Count++;
        }
        return Count;
    }


    /** Divide inputNum by digits , Its will work only for integer number!!! */
    template<typename INT_T>
    void DivideNumsByDig(UINT64 InputNum, INT_T* (&Array))
    {
        UINT32 len = CountDigits(InputNum);

        for(UINT32 i=0; i<len; i++)
        {
            if(i==0)
            {
                Array[i] = InputNum/pow(10,len-1);
                continue;
            }
            Array[i] = InputNum/pow(10,(len-1)-i);
            Array[i] %= 10;
        }
    }

    /** Divide inputNum(type string) by digits numbers and convert to integer type */
    template<typename INT_T>
    void DivideNumsByDig(std::string InputNum, INT_T* (&Array))
    {
        UINT32 len = InputNum.length();
        for(UINT32 i=0; i<len; i++)
        {
            CHAR temp = InputNum[i];
            Array[i] = atoi(&temp); // char to int
        }
    }

    /** Get maximum of number in the array */
    template <typename Type>
    Type GetMax(const Type* Array, UINT32 Size)
    {
        Type Max = Array[0];

        for(UINT32 i=0; i<Size; i++)
        {
            if(Max<Array[i])
                Max = Array[i];
        }
        return Max;
    }

    /** Get maximum of number in the range: StartIterator...EndIterator */
    template <typename FowardIterator>
    FowardIterator GetMax(FowardIterator* Start, FowardIterator* End)
    {
        FowardIterator* Max = Start;
        UINT32 Size = End - Start;

        for(UINT32 i=0; i<=Size; i++)
        {
            if(*(Start+i) > *Max)
            *Max = *(Start+i);
        }
        return *Max;
    }

    /** Get maximum of 2 numbers */
    template <typename Type>
    inline Type GetMax(Type Num1, Type Num2)
    {
        return (Num1>=Num2 ? Num1:Num2);
    }

    /** Get minimum of number in the array */
    template <typename Type>
    Type GetMin(const Type* Array, UINT32 Size)
    {
        Type Min = Array[0];

        for(UINT32 i=0; i<Size; i++)
        {
            if(Min>Array[i])
                Min = Array[i];
        }
        return Min;
    }

    /** Get minimum of number in the range: StartIterator...EndIterator */
    template <typename FowardIterator>
    FowardIterator GetMin(FowardIterator* Start, FowardIterator* End)
    {
        FowardIterator* Min = Start;
        UINT32 Size = End - Start;

        for(UINT32 i=0; i<=Size; i++)
        {
            if(*(Start+i) < *Min)
                *Min = *(Start+i);
        }
        return *Min;
    }

    /** Get minimum of 2 numbers */
    template <typename Type>
    inline Type GetMin(Type Num1, Type Num2)
    {
        return (Num1<=Num2 ? Num1:Num2);
    }


    /** Sort 1D array by increasing order(Algorithm Bubble Sort, complex O(N^2)) */
    template <typename ArrayType>
    void SortInc_B(ArrayType* (&Array), UINT32 Size)
    {
        ArrayType Temp; // temporary variable

        //Algorithm "Buble sort"
        for(UINT32 i=0; i<Size; i++)
        {
            for(UINT32 j=0; j<(Size-i-1); j++)
            {
                if(Array[j]>Array[j+1])
                {
                    Temp = Array[j];
                    Array[j] = Array[j+1];
                    Array[j+1]= Temp;
                }
            }
        }
    }

    /** Sort 1D array by decreasing order(Algorithm Bubble Sort, complex O(N^2)) */
    template <typename ArrayType>
    void SortDec_B(ArrayType *(&Array), UINT32 Size) // Sort array by decreasing order
    {
        ArrayType Temp; // temporary variable

        //Algorithm "Buble sort"
        for(UINT32 i=0; i<Size; i++)
        {
            for(UINT32 j=0; j<(Size-i-1); j++)
            {
                if(Array[j]<Array[j+1])
                {
                    Temp = Array[j];
                    Array[j] = Array[j+1];
                    Array[j+1]= Temp;
                }
            }
        }
    }


    /** Input whatever array in the Console */
    template <typename Type>
    inline void ScanArray(Type* Array, UINT32 Size)
    {
        for(UINT32 i=0; i<Size; i++)
        {
            std::cin >> Array[i];
        }
    }

    /** Print whatever array to the display */
    template <typename Type>
    inline void PrintArray(const Type* Array, UINT32 Size)
    {
        for(UINT32 i=0; i<Size; i++)
        {
            std::cout << Array[i] << std::setw(4);
        }
    }

    /** Recursive Binary Search complex O(log(n)) */
    template <typename Type>
    Type Recursive_BinarySearch(const Type* Array, UINT32 Left, UINT32 Right, const Type FindElement)
    {
        //INT32 Left = 1;
        //INT32 Right = Size;
        UINT32 Middle = (Left + Right)/2;
        bool Founded = false;

        if(Left>Right)
        {
            return -1;
        }
        else if(Left<=Right)
        {
            Middle = (Left + Right)/2;

            if(Array[Middle]==FindElement)
            {
                return Middle;
                Founded = true;
            }
            else if(Array[Middle]>FindElement)
            {
                return Recursive_BinarySearch(Array, Left, Right-1, FindElement);
            }
            else if(Array[Middle]<FindElement)
            {
                return Recursive_BinarySearch(Array, Left, Right+1, FindElement);
            }
        }
    }

    /** Simple Binary Search complex O(log(n)), if element was found in the sorted array, then function returned position of element in sorted array, otherwise function will return -1 */
    template <typename Type>
    Type BinarySearch(const Type* Array, UINT32 Size, Type FindElement)
    {
        INT32 low = 1;
        INT32 high = Size;
        bool Founded = false;

        while(low<=high)
        {
            INT32 mid = (low+high)/2;

            if(Array[mid]>FindElement)
            {
                high = mid-1;
            }
            else if(Array[mid]<FindElement)
            {
                low = mid+1;
            }
            else if(Array[mid]==FindElement)
            {
                return mid;
                Founded = true;
            }
        }
        if(!Founded)
        {
            return -1;
        }
    }

    /** Function calculate the greatest common divisor(GCD) two natural numbers */
    template<typename INT_T>
    INT_T GCD(INT_T num1, INT_T num2)
    {
        while(num1*num2>0)
        {
            if(num1>=num2)
            {
                num1 %= num2;
            }
            else{
                num2 %= num1;
            }
        }
        return num1+num2; //���
    }

    /** Function calculate the least common multiples(LCM) two natural numbers */
    template<typename INT_T>
    INT_T LCM(INT_T num1, INT_T num2)
    {
        return (num1*num2)/GCD(num1, num2); //��� (���������� ����� �������)
    }

//-------------------------------------------------------------------------------------------------------------------------------
//----------------------------------------------  CLASSES  ----------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------------------------
	struct HashEntry
	{
		typedef std::string String;
		String Value;

		HashEntry(String Value)
		{
			this->Value = Value;
		}
	};
//-------------------------------------------------------------------------------------------------------------------------------
    /** class Hash Table */
	class HashTable
    {
    private:
		const INT16 TABLE_SIZE = 127; // ������������ ������ �����
		typedef std::string String;
		HashEntry **Table; //������� �����

		INT64 HashFunction(String Value)
		{
			INT64 Sum = 0; //ASCII sum
			for (auto i : Value)
			{
				Sum += i;
			}
			return Sum % TABLE_SIZE;
		}
    public:
		HashTable()
		{
			Table = new HashEntry*[TABLE_SIZE];
			for (INT16 i = 0; i < TABLE_SIZE; i++)
			{
				Table[i] = nullptr;
			}
		}
		HashTable(const HashTable& OtherHashTable)
		{
			//this->TABLE_SIZE = OtherHashTable.TABLE_SIZE;
			this->Table = OtherHashTable.Table;
			for (INT16 i = 0; i < this->TABLE_SIZE; i++)
			{
				this->Table[i] = OtherHashTable.Table[i];
			}
		}
        ~HashTable()
        {
			for (INT16 i = 0; i < TABLE_SIZE; i++)
			{
				delete Table[i];
			}
			delete Table;
        }

		void insert(String Value)
		{
			INT64 HashNumber = HashFunction(Value);
			if (Table[HashNumber] == nullptr)
			{
				Table[HashNumber] = new HashEntry(Value);
			}
		}

		void remove(String Value)
		{
			INT64 HashNumber = HashFunction(Value);
			if (Table[HashNumber] != nullptr)
			{
				Table[HashNumber] = nullptr;
			}
		}

		HashEntry* find(String Value)
		{
			INT64 HashNumber = HashFunction(Value);
			if (Table[HashNumber] == nullptr)
			{
				return nullptr;
			}
			else{
				return Table[HashNumber];
			}
		}

		INT64 getHash(String Value)
		{
			INT64 HashNumber = HashFunction(Value);
			return HashNumber;
		}
    };

//-------------------------------------------------------------------------------------------------------------------------------

    /** Class Point(Space and Flat) */
    template<class Type>
    class Point
    {
    private:
        //Coordinates of Point in the Space
        Type cordX;
        Type cordY;
        Type cordZ;

    public:
        //Constructors
        Point(): cordX(0), cordY(0), cordZ(0){}
        Point(Type X, Type Y, Type Z=0): cordX(X), cordY(Y), cordZ(Z){}
        ~Point(){} //Destructor

        void ShowCoords()
        {
            Type X = cordX;
            Type Y = cordY;
            Type Z = cordZ;
            std::cout <<"X="<< X <<" Y="<< Y <<" Z="<< Z <<"\n";
        }
    };

//-------------------------------------------------------------------------------------------------------------------------------

    template <class T>
    class Determinant
    {
    private:
        T det; // element of determinant
        T result;
        UINT32 size_det;

        void GetMinor(T** (&Array), T** p, UINT16 i, UINT16 j, UINT16 m)
        {
            UINT16 ki=0, kj=0, di=0, dj=0;

            for (; ki<(m-1); ki++) // control index of row
            {
                if (ki==i) di=1;
                dj=0;
                for (; kj<(m-1); kj++) // control index of col
                {
                    if (kj == j) dj=1;
                    p[ki][kj] = det[ki+di][kj+dj];
                }
            }
        }

    public:
		Determinant() {} // constructor
		~Determinant() {} // destructor

        void Display(T** (&Array), UINT32 Size)
        {
            //**                    **//
                using namespace std;
            //**                    **//

            Size = size_det;
            Array = det;

            cout << "\n\n" << "|";
            for(UINT32 i=0; i<size_det; i++)
            {
                UINT32 j=0;
                for(; j<size_det; j++)
                {
                    cout << setw(4) << det[i][j] << setw(4);
                }
                cout << "|" << endl;

                if((i != size_det-1) && (j != size_det-1)){
                    cout << "|";
                }
            }
            cout << "\n\n";
        }

        T Compute(T** (&Array), UINT32 Size) // Recruisive computing determinant
        {
            Size = size_det;
            Array = det;

            T **p = nullptr;

            INT16 k = 1; //(-1)^(i+j)
            UINT32 i=0,j=0;

            p = new T*[Size];

            for (; i<Size; i++){
                p[i] = new T[Size];
            }

            if (size_det == 1)
            {
                result = det[0][0];
                return result;
            }
            else if(size_det == 2)
            {
                result = (det[0][0]*det[1][2]) - (det[0][1]*det[1][0]);
                return result;
            }
            else if(size_det > 2)
            {
                for (; i<size_det; i++)
                {
                    GetMinor(det, p, i, 0, size_det);
                    //cout << det[i][j] << endl;
                    result = result + (k*(det[i][0])*ComputeDeterminant(p,(size_det-1)));
                    k = -k;
                }
            }
        }

        void ShowResult(T** (&array), UINT32 Size)
        {
            std::cout << "D=" << result;
        }
    };
//-------------------------------------------------------------------------------------------------------------------------------
    /** class Big Integer */
    class BigInteger
    {
    private:
        UINT16* ArrayDigits;
        UINT16 SIZE;
        bool isCorrectNumber = true;

    public:
        const static UINT16 ONE = 1;
        const static UINT16 ZERO = 0;

        BigInteger()
        {
            SIZE = 1;
            ArrayDigits = new UINT16[SIZE];
            ArrayDigits[0] = this->ZERO;

        }

        BigInteger(const BigInteger *A)
        {
            SIZE = A->SIZE;
            ArrayDigits = new UINT16[SIZE];
            for(UINT64 i=0; i<SIZE; i++)
            {
                ArrayDigits[i] = A->ArrayDigits[i];
            }
        }

        BigInteger(std::string BigIntegerNumber)
        {
            SIZE = BigIntegerNumber.length();
            ArrayDigits = new UINT16[SIZE];

            for(INT16 i=SIZE-1; i>=0; i--)
            {
                CHAR temp = BigIntegerNumber[(SIZE-1)-i];

                try
                {
                    if(!isdigit(temp))
                    {
                        isCorrectNumber = false;
                        throw std::runtime_error("\nConvert Error Big Integer \n");
                    }
                    else{
                        UINT16 num = atoi(&temp); // char to int

                        if(num>=10)
                        {
                            num /= 10;
                        }
                        ArrayDigits[i] = num;
                    }
                }
                catch(std::runtime_error& e)
                {
                    std::cerr << e.what();
                    break;
                }
            }
        }

        ~BigInteger()
        {
            if(ArrayDigits != nullptr)
            {
                std::cout<<"Destructor called\n";
                delete[] ArrayDigits;
            }
        }

        std::string ToString()
        {
            std::string str;
            std::string temp;

            for(INT64 i=SIZE-1; i>=0 && isCorrectNumber; i--)
            {
                temp = std::to_string(ArrayDigits[i]);
                //str += temp;
                str.push_back(temp[0]);
                std::cout << str << "\n";
            }
            return str;
        }

        UINT64 getSize()
        {
            return SIZE;
        }

        BigInteger& operator=(const BigInteger*)
        {
            return *this;
        }

        BigInteger& operator=(std::string BigIntegerNumber)
        {
            SIZE = BigIntegerNumber.length();
            ArrayDigits = new UINT16[SIZE];

            for(INT16 i=SIZE-1; i>=0; i--)
            {
                CHAR temp = BigIntegerNumber[(SIZE-1)-i];

                try
                {
                    if(!isdigit(temp))
                    {
                        isCorrectNumber = false;
                        throw std::runtime_error("\nConvert Error Big Integer \n");
                    }
                    else{
                        UINT16 num = atoi(&temp); // char to int

                        if(num>=10)
                        {
                            num /= 10;
                        }
                        ArrayDigits[i] = num;
                    }
                }
                catch(std::runtime_error& e)
                {
                    std::cerr << e.what();
                    break;
                }
            }
            return *this;
        }

        BigInteger& operator+(const BigInteger &A)
        {
            BigInteger* res = new BigInteger(); //res - result
            UINT16 carry = 0;

            if(SIZE == A.SIZE)
            {
                for(UINT64 i=0; i<SIZE; i++)
                {
                    UINT16 sum = ArrayDigits[i] + A.ArrayDigits[i] + carry;
                    if(sum>=10)
                    {
                        res->ArrayDigits[i] = sum - 10;
                        carry = 1;
                    }
                    else{
                        res->ArrayDigits[i] = sum;
                        carry = 0;
                    }
                }
                res->SIZE = SIZE;
            }
            return *res;
        }




        #if DEBUG >= 1
        void test_ArrayDigits()
        {
            for(UINT64 i=0; i<SIZE; i++)
            {
                std::cout<<ArrayDigits[i]<<"\n";
            }
        }
        #endif // DEBUG

    };

//-------------------------------------------------------------------------------------------------------------------------------
    /** Template of Class Stack, Stack is structure data that realized the interface of "First-In-Last-Out" (FILO)  */
    template<class Type>
    class Stack
    {
    private:
        UINT64 Size; //������ �����
        UINT64 Count; //���������� ����������� ��������� � �����
        Type* Array; //������� ������
        Type* NewArray = nullptr; //��������������� ������

    public:
        Stack(UINT64 MaxSize = 0) //����������� �� ���������
        {
            Size = MaxSize;
            Count = 0;
            if(MaxSize==0)
                Array = new Type[++MaxSize]{0};
            else
                Array = new Type[MaxSize]{0};
        }
        Stack(const Stack<Type> &OtherStack) //���������� �����������
        {
            Size = OtherStack.Size;
            Count = OtherStack.Count;
            Array = new Type[Size];

            for(UINT64 i=0; i<Size; i++)
            {
                Array[i] = OtherStack.Array[i];
            }
        }
        ~Stack() //����������
        {
            if(Array != nullptr)
                delete[] Array;
            if(NewArray != nullptr)
                delete[] NewArray;
        }
        void push_back(const Type Value) //��������� �������� � ����� �����
        {
            NewArray = new Type[++Size]{0};
            for(UINT64 i=0; i<Size; i++)
            {
                NewArray[i] = Array[i];
            }
            delete[] Array;
            Array = NewArray;

            Array[Count++] = Value;
        }
        void pop_back() //������� �������� �� ����� �����
        {
            try
            {
                if(Size==0 || Count==0)
                    throw std::runtime_error("\nERROR: the Size of Stack is NULL \n");

                NewArray = new Type[Size]{0};
                for(UINT64 i=0; i<Size; i++)
                {
                    NewArray[i] = Array[i];
                }
                delete[] Array;
                Array = NewArray;

                Size--;
                Count--;
                Array[Count] = 0;
            }
            catch(std::runtime_error& e)
            {
                std::cerr << e.what();
            }
        }
        Type back() //�������� ������ � ����� �����
        {
            if(Count!=0)
                return Array[Count-1];
            else
                return 0;
        }
        UINT64 size() const //������ ����������� ��������� � �����
        {
            return Count;
        }
        /*UINT64 capacity() const //����������� ������ �����
        {
            return Size;
        }
        */
        bool empty() //��������� �����
        {
            if(Size>=1 && Count>0)
                return false;
            else
                return true;
        }
        void clear() //������� ��� �������� �� �����
        {
            for(UINT64 i=0; i<Size; i++)
            {
                Array[i] = 0;
            }
            Count = 0;
            Size = 0;
        }

        Stack<Type>& operator=(const Stack<Type>& OtherStack) const
        {
            this->Size = OtherStack.Size;
            this->Count = OtherStack.Count;
            Array = new Type[Size];

            for(UINT64 i=0; i<Size; i++)
            {
                this->Array[i] = OtherStack.Array[i];
            }
            return *this;
        }

        #if DEBUG >= 1
        void print()
        {
            for(UINT64 i=0; i<Size; i++)
            {
                std::cout << Array[i] <<" ";
            }
        }
        #endif // DEBUG
    };

//-------------------------------------------------------------------------------------------------------------------------------
    /** Template of Class Queue, Queue is structure data that realized the interface of "First-In-First-Out" (FIFO)  */
    template<class Type>
    class Queue
    {
    private:
        UINT64 Size;    //������ �������
        UINT64 Begin;   //������ �������
        UINT64 End;     //����� �������
        Type* Array;    //������� ������
        Type* NewArray = nullptr; //��������������� ������

    public:
        Queue(UINT64 MaxSize = 0) //����������� �� ���������
        {
            Size = MaxSize;
            Begin = 0;
            End = 0;
            if(MaxSize==0)
                Array = new Type[++MaxSize]{0};
            else
                Array = new Type[MaxSize]{0};
        }
        Queue(const Queue<Type> &OtherQueue) //���������� �����������
        {
            Size = OtherQueue.Size;
            Begin = OtherQueue.Begin;
            End = OtherQueue.End;
            Array = new Type[Size];

            for(UINT64 i=0; i<Size; i++)
            {
                Array[i] = OtherQueue.Array[i];
            }
        }
        ~Queue() //����������
        {
            if(Array != nullptr)
                delete[] Array;
            if(NewArray != nullptr)
                delete[] NewArray;
        }
        void push_back(const Type Value) //��������� �������� � ����� �������
        {
            NewArray = new Type[++Size]{0};
            for(UINT64 i=0; i<Size; i++)
            {
                NewArray[i] = Array[i];
            }
            delete[] Array;
            Array = NewArray;

            End = Size - 1;
            Array[End] = Value;
        }
        void pop_front() //��������� �������� �� �������
        {
            try
            {
                if(Size==0 || End<0 || Begin>End)
                    throw std::runtime_error("\nERROR: the Size of Queue is NULL \n");

                NewArray = new Type[Size]{0};
                for(UINT64 i=0; i<Size; i++)
                {
                    if(i==Size-1)
                    {
                        NewArray[i] = 0; //������� �������� ������� � �������
                    }
                    else{
                        NewArray[i] = Array[i+1]; //����� ��������� � �����(<--)
                    }
                }

                delete[] Array;
                Array = NewArray;

                Size--;
                End = Size - 1;
                //Begin = 0;
            }
            catch(std::runtime_error& e)
            {
                std::cerr << e.what();
            }

        }
        Type front() //�������� ������ � ������ �������
        {
            if(End>=Begin)
                return Array[Begin];
            else
                return 0;
        }
        Type back() //�������� ������ � ����� �������
        {
            if(End>=Begin)
                return Array[End];
            else
                return 0;
        }
        UINT64 size() const //������ ����������� ��������� � �������
        {
            return Size;
        }
        bool empty() //��������� �������
        {
            if(Size>=1 && End>=Begin)
                return false;
            else
                return true;
        }
        void clear() //������� ��� �������� �� �������
        {
            for(UINT64 i=0; i<Size; i++)
            {
                Array[i] = 0;
            }
            Begin = 0;
            End = 0;
            Size = 0;
        }

        Queue<Type>& operator=(const Queue<Type>& OtherQueue)
        {
            this->Size = OtherQueue.Size;
            this->Begin = OtherQueue.Begin;
            this->End = OtherQueue.End;
            Array = new Type[Size];

            for(UINT64 i=0; i<Size; i++)
            {
                this->Array[i] = OtherQueue.Array[i];
            }
            return *this;
        }

        #if DEBUG >= 1
        void print()
        {
            std::cout << "Begin: "<<Begin<<" End: "<<End<<"\n";
            for(UINT64 i=0; i<Size; i++)
            {
                std::cout << Array[i] <<" ";
            }
            //std::cout <<"\n"<< Array[End+1] <<" \n";
        }
        #endif // DEBUG
    };

//-------------------------------------------------------------------------------------------------------------------------------
    /** Template of Class Deque, Deque is structure data that realized the interface of "Double-Ended-Queue"  */
    template<class Type>
    class Deque
    {
    private:
        UINT64 Size;                //������ ����
        UINT64 Begin;               //������ ���� (Iterator begin)
        UINT64 End;                 //����� ���� (Iterator end)
        //UINT64 IndexOfElement;      //������ ��������
        Type* Array;                //������� ������
        Type* NewArray = nullptr;   //��������������� ������

        typedef Type* iterator;
        typedef Type& reference;
        typedef const Type* const_iterator;
        typedef const Type& const_reference;

    public:
        Deque(UINT64 MaxSize = 0) //����������� �� ���������
        {
            Size = MaxSize;
            Begin = 0;
            End = 0;
            //IndexOfElement = 0;
            if(MaxSize==0)
                Array = new Type[++MaxSize]{0};
            else
                Array = new Type[MaxSize]{0};
        }
        Deque(const Deque<Type> &OtherDeque) //���������� �����������
        {
            Size = OtherDeque.Size;
            Begin = OtherDeque.Begin;
            End = OtherDeque.End;
            Array = new Type[Size];

            for(UINT64 i=0; i<Size; i++)
            {
                Array[i] = OtherDeque.Array[i];
            }
        }
        ~Deque() //����������
        {
            if(Array != nullptr)
                delete[] Array;
            if(NewArray != nullptr)
                delete[] NewArray;
        }
        void push_back(const Type Value) //��������� �������� � ����� ����
        {
            NewArray = new Type[++Size]{0};
            for(UINT64 i=0; i<Size; i++)
            {
                NewArray[i] = Array[i];
            }
            delete[] Array;
            Array = NewArray;

            End = Size - 1;
            Array[End] = Value;
        }
        void push_front(const Type Value) //��������� �������� � ������ ����
        {
            NewArray = new Type[++Size]{0};
            for(UINT64 i=0; i<Size; i++)
            {
                if(i!=Size-1)
                {
                    NewArray[i+1] = Array[i]; //����� ��������� ������� � ������(-->)
                }
            }
            delete[] Array;
            Array = NewArray;

            End = Size - 1;
            Array[Begin] = Value;   //Begin = 0
                                    //��������� �������� � ������ �������
        }
        void pop_front() //��������� �������� �� ������ ����
        {
            try
            {
                if(Size==0 || End<0 || Begin>End)
                    throw std::runtime_error("\nERROR: the Size of Queue is NULL \n");

                NewArray = new Type[Size]{0};
                for(UINT64 i=0; i<Size; i++)
                {
                    if(i==Size-1)
                    {
                        NewArray[i] = 0; //������� �������� ������� � �������
                    }
                    else{
                        NewArray[i] = Array[i+1]; //����� ��������� ������� � �����(<--)
                    }
                }

                delete[] Array;
                Array = NewArray;

                Size--;
                End = Size - 1;
                //Begin = 0;
            }
            catch(std::runtime_error& e)
            {
                std::cerr << e.what();
            }
        }

        void assign(const Type Value, const UINT_T Pos)
        {
            try
            {
                if(Pos>=Size)
                {
                    throw std::runtime_error("\nERROR: the index of element is largest than Deque size \n");
                }
                this->Array[Pos] = Value;
            }
            catch(std::runtime_error& e)
            {
                std::cerr << e.what();
                return;
            }
        }

        void pop_back() //������� �������� �� ����� ����
        {
            try
            {
                if(Size==0 || End<0 || Begin>End)
                    throw std::runtime_error("\nERROR: the Size of Deque is NULL \n");

                NewArray = new Type[Size]{0};
                for(UINT64 i=0; i<Size; i++)
                {
                    NewArray[i] = Array[i];
                }
                delete[] Array;
                Array = NewArray;

                Size--;
                End = Size - 1;
                //Array[Count] = 0;
            }
            catch(std::runtime_error& e)
            {
                std::cerr << e.what();
            }
        }
        const_reference front() const //�������� ������ � ������ ����
        {
            //if(End>=Begin)
                return *(this->Array+this->Begin);
        }
        const_reference back() const //�������� ������ � ����� ����
        {
            //if(End>=Begin)
                return *(this->Array+this->End);
        }
        UINT64 size() const //������ ����������� ��������� � ����
        {
            return Size;
        }
        bool empty() //��������� ����
        {
            if(Size>=1 && End>=Begin)
                return false;
            else
                return true;
        }

        void resize(UINT_T NewSize)
        {
            try
            {
                if(NewSize<=0)
                    throw std::runtime_error("\nERROR: the Size of Deque is least than ZERO \n");

                NewArray = new Type[NewSize]{0};
                if(NewSize>=Size)
                {
                    for(UINT64 i=0; i<Size; i++)
                    {
                        NewArray[i] = Array[i];
                    }
                }
                else{
                    for(UINT64 i=0; i<NewSize; i++)
                    {
                        NewArray[i] = Array[i];
                    }
                }
                delete[] Array;
                Array = NewArray;

                Size = NewSize;
                End = Size - 1;

            }
            catch(std::runtime_error& e)
            {
                std::cerr << e.what();
            }
        }

        void clear() //������� ��� �������� �� ����
        {
            for(UINT64 i=0; i<Size; i++)
            {
                Array[i] = 0;
            }
            Begin = 0;
            End = 0;
            Size = 0;
            //IndexOfElement = 0;
        }

        const_iterator begin() const // �������� �� ������ ����
        {
            return this->Array + this->Begin;
        }
        const_iterator end() const // �������� �� ����� ����
        {
            return this->Array + this->End;
        }

        const_reference at(UINT_T ElementIndex) const
        {
            try
            {
                //IndexOfElement = ElementIndex;
                if(ElementIndex>=Size)
                {
                    throw std::runtime_error("\nERROR: the index of element is largest than Deque size \n");
                }
                return *(this->Array+ElementIndex);
            }
            catch(std::runtime_error& e)
            {
                std::cerr << e.what();
                return 0;
            }
        }

        const_reference operator[](UINT_T ElementIndex) const
        {
            //IndexOfElement = ElementIndex;
            return *(this->Array+ElementIndex);
        }

        Deque<Type>& operator=(const Deque<Type>& OtherDeque)
        {
            this->Size = OtherDeque.Size;
            this->Begin = OtherDeque.Begin;
            this->End = OtherDeque.End;

            Array = new Type[this->Size];
            for(UINT64 i=0; i<this->Size; i++)
            {
                Array[i] = OtherDeque.Array[i];
            }

            return *this;
        }

        /*const void operator=(const Type Value) const
        {
            //this->Array[IndexOfElement] = Value;
            std::cout << "void operator= \n";
            return;
        }
        */


        #if DEBUG >= 1
        void print()
        {
            std::cout << "Begin: "<<Begin<<" End: "<<End<<"\n";
            for(UINT64 i=0; i<Size; i++)
            {
                std::cout << Array[i] <<" ";
            }
            //std::cout <<"\n"<< Array[End+1] <<" \n";
        }
        #endif // DEBUG
    };
//-------------------------------------------------------------------------------------------------------------------------------
    /** Template struct of Node, Node is structure data which saved itself data and address of next and previous data */
    template<class Type>
    struct Node
    {
        Type Data;          //data
        Node<Type>* Next;   //pointer to next data
        Node<Type>* Prev;   //pointer to previous data

        Node<Type>* Right;  //pointer to right node (used for binary tree)
        Node<Type>* Left;   //pointer to left node (used for binary tree)

        Node(): Next(nullptr), Prev(nullptr), Right(nullptr), Left(nullptr) {} //Default constructor
        Node(const Type& setData): Data(setData), Next(nullptr), Prev(nullptr), Right(nullptr), Left(nullptr) {}
        Node(const Node<Type>& OtherNode) //Copying constructor
        {
            this->Data = OtherNode.Data;
            this->Left = OtherNode.Left;
            this->Right = OtherNode.Right;
            this->Next = OtherNode.Next;
            this->Prev = OtherNode.Prev;
        }
    };
//-------------------------------------------------------------------------------------------------------------------------------
    /** Template class of List                                                             */
    /** Forward List is structure data which saved itself data and address of next data    */
    /**                                                                                    */
    template<class Type>
    class ForwardList
    {
    private:
        Node<Type>* Head;     //��������� �� ������ ������
        typedef Node<Type>* iterator;
        typedef const Node<Type>* const_iterator;
        typedef Type& reference;
        typedef const Type& const_reference;

    public:
        ForwardList(): Head(nullptr) {}
        ForwardList(const ForwardList<Type>& OtherForwardList)
        {
            this->Head = OtherForwardList.Head;
        }
        ~ForwardList() //����������
        {
            while(Head!=nullptr)  //���� �� ������ �� �����
            {
                Node<Type>* Temp = Head->Next;  //��������� ���������� ��� �������� ������ ���������� ��������
                delete Head;                    //����������� ����� ������������ ������
                Head = Temp;                    //������ ����� �� ���������
            }
        }

        void push_front(const Type Value) //������� ���������� ��������� � ������
        {
            Node<Type>* Temp = new Node<Type>(); //��� ������ ������ ���������� ������
            Temp->Data = Value; //���������� Valu � ������� ���������  Node
            Temp->Next = Head;  //���������, ��� ����. ������� ��� ������ �� ������ Head
            Head = Temp;        //���������, ��� ��������� �������� ������� ��� ������ ��� ���������
        }

        void push_back(const Type Value) //������� ���������� ��������� � ������
        {
            Node<Type>* Temp = new Node<Type>(); //��� ������ ������ ���������� ������
            Node<Type>* NewNode = new Node<Type>(Value);
            Temp = Head;
            while(Temp->Next != nullptr)
            {
                Temp = Temp->Next;
            }
            Temp->Next = NewNode;
        }

        void assign(const Type Value, const UINT_T Pos)
        {
            if(Head!=nullptr && Pos <= this->size()) //������ �������� �� �� ��� ������ �� ���� � Pos �� ��������� ����� ��� ���������
            {
                Node<Type>* Temp = Head;
                for(UINT64 i=1; i<=Pos; i++)
                {
                    Temp = Temp->Next; //������ ����� N ���
                }
                Temp->Data = Value;
            }
        }

        void insert(const Type Value, const UINT_T Pos)
        {
            if(Head!=nullptr && Pos!=1 && Pos < this->size() ) //������ �������� �� �� ��� ������ �� ���� � Pos �� ��������� ����� ��� ���������
            {
                Node<Type>* BeforeNode = Head;
                Node<Type>* AfterNode = Head;
                Node<Type>* NewNode = new Node<Type>(Value);

                for(UINT64 i=1; i<Pos-1; i++)
                {
                    BeforeNode = BeforeNode->Next;
                }
                for(UINT64 i=1; i<Pos; i++)
                {
                    AfterNode = AfterNode->Next;
                }
                BeforeNode->Next = NewNode;
                NewNode->Next = AfterNode;

                //std::cout << BeforeNode <<"\n";
                //std::cout << AfterNode <<"\n";
                //if(BeforeNode==AfterNode) puts("YES");
            }
            else if(Pos==1 && this->size()>=1)
            {
                this->push_front(Value);
            }
            else if(Pos==this->size())
            {
                this->push_back(Value);
            }
        }

        void erase(const Type Element) //�������� Element �� ������
        {
            Node<Type>* ElementNode = Head;
            Node<Type>* Temp = Head;

            //����� Element �� ������
            while(ElementNode->Next!=nullptr && ElementNode->Data != Element)
            {
                ElementNode = ElementNode->Next;
            }

            //������� ElementNode �� ������
            if(Temp==ElementNode)
            {
                Temp->Next = ElementNode->Next;
            }
            else{
                while(Temp->Next != ElementNode)
                {
                    Temp = Temp->Next;
                }
                Temp->Next = ElementNode->Next;
            }
        }

        void pop_front() //�������� ������ ������ �� ������ ������
        {
            Head = Head->Next;
        }

        void pop_front_nodes(const UINT_T N) //�������� N ���� �� ������ ������
        {
            if(Head!=nullptr  &&  N <= this->size()) //������ �������� �� �� ��� ������ �� ���� � N �� ��������� ����� ��� ���������
            {
                for(UINT64 i=1; i<=N; i++)
                {
                    this->pop_front(); //������ ����� N ���
                }
            }
        }

        void pop_back() //�������� ������ ������ �� ����� ������
        {
            Node<Type>* Temp = Head;
            Node<Type>* LastElement = nullptr;

            while(Temp->Next!=nullptr) //���� � ����� ������
            {
                Temp = Temp->Next;
            }
            LastElement = Temp;
            LastElement->Data = Temp->Data; //�������� ������ �� ���������� ������ � ����� LastElement
            this->erase(LastElement->Data); //����� �������� ���������� �������� �� ������
        }

        void pop_back_nodes(const UINT_T N) //�������� N ���� �� ����� ������
        {
            if(Head!=nullptr  &&  N <= this->size()) //������ �������� �� �� ��� ������ �� ���� � N �� ��������� ����� ��� ���������
            {
                for(UINT64 i=1; i<=N; i++)
                {
                    this->pop_back();
                }
            }
        }

        const_reference front() const
        {
            return Head->Data;
        }

        const_reference back() const
        {
            Node<Type>* Temp = Head;
            while(Temp->Next!=nullptr)
            {
                Temp = Temp->Next;
            }
            return Temp->Data;
        }

        const_iterator begin() const
        {
            return Head;
        }

        const_iterator end() const
        {
            Node<Type>* Temp = Head;
            while(Temp->Next!=nullptr)
            {
                Temp = Temp->Next;
            }
            return Temp;
        }

        ForwardList<Type>& operator=(const ForwardList<Type>& OtherForwardList) const
        {
            this->Head = OtherForwardList.Head;
            return *this;
        }

        void clear()
        {
            while (Head!=nullptr)  //���� �� ������ �� �����
            {
                Node<Type>* Temp = Head->Next;  //��������� ���������� ��� �������� ������ ���������� ��������
                delete Head;                    //����������� ����� ������������ ������
                Head = Temp;                    //������ ����� �� ���������
            }
        }

        bool empty()
        {
            if(this->size())
                return false;
            else
                return true;
        }

        UINT64 size()
        {
            UINT64 Count = 0;
            Node<Type>* Temp = Head;
            while(Temp!=nullptr)
            {
                Count++;
                Temp = Temp->Next;
            }
            return Count;
        }

        #if DEBUG >= 1
        void print()
        {
            Node<Type>* Temp = Head; //���������� ���������, ������� ���������� �� ����� ������ ������ ������
            while(Temp!=nullptr) //�� ��� ��� ���� �� �������� ������ ��������
            {
                std::cout << Temp->Data <<" -> ";
                Temp = Temp->Next; //���������, ��� ����� ��� ����� ��������� �������
            }
        }
        #endif // DEBUG

    };
//-------------------------------------------------------------------------------------------------------------------------------
    /** Template class of Binary tree, Binary Tree is data structure, which trees are hierarchical  */
    /** They consist of interconnected nodes, The topmost node is called the root.                  */
    /** This Unique Binary Tree                                                                     */
    template <class Type>
    class BinaryTree
    {
    private:
        Node<Type>* Root;

    public:
        BinaryTree(): Root(nullptr) {}
        BinaryTree(const Type Value)
        {
            Root->Data = Value;
            Root = nullptr;
        }
        BinaryTree(const BinaryTree<Type>& OtherBinaryTree)
        {
            this->Root = OtherBinaryTree.Root;
        }
        ~BinaryTree()
        {
            if(this->Root!=nullptr)
            {
                this->clear(this->getRootTree());
            }
        }

        void insert(const Type Item)
        {
            Node<Type>* Temp1 = new Node<Type>();
            if(Root==nullptr)
            {
                Root = Temp1;
                Root->Data = Item;
                return;
            }

            Temp1 = Root;
            Node<Type>* Temp2 = nullptr;

            while(Temp1 != nullptr) //������ �� ������� �������
            {
                if(Temp1->Data == Item)
                {
                    //Temp1->Data = Item;
                    return;
                }
                else{
                    Temp2 = Temp1;
                    if(Item < Temp1->Data)
                    {
                        Temp1 = Temp1->Left;
                    }
                    else{
                        Temp1 = Temp1->Right;
                    }
                }
            }

            Node<Type>* NewNode = new Node<Type>(Item);
            if(Temp2 == nullptr)
            {
                Root = Temp2;
            }
            else{
                if(Item < Temp2->Data)
                {
                    Temp2->Left = NewNode;
                }
                else{
                    Temp2->Right = NewNode;
                }
            }
        }

        void remove(const Type Item)
        {
            Node<Type>* Temp1 = Root;
            Node<Type>* Temp2 = nullptr;

            while(Temp1 != nullptr)
            {
                if(Temp1->Data == Item)
                {
                    //Temp1->Data = Item;
                    break;
                }
                else{
                    Temp2 = Temp1;
                    if(Item < Temp1->Data)
                    {
                        Temp1 = Temp1->Left;
                    }
                    else{
                        Temp1 = Temp1->Right;
                    }
                }
            }

            if(Temp1==nullptr)
            {
                return;
            }

            if(Temp1->Right==nullptr)
            {
                if(Temp2==nullptr)
                {
                    Root = Temp1->Left;
                }
                else{
                    if(Temp1 == Temp2->Left)
                    {
                        Temp2->Left = Temp1->Left;
                    }
                    else{
                        Temp2->Right = Temp1->Left;
                    }
                }
            }
            else{
                Node<Type>* LeftMost = Temp1->Right;
                Temp2 = nullptr;

                while(LeftMost->Left!=nullptr)
                {
                    Temp2 = LeftMost;
                    LeftMost = LeftMost->Left;
                }

                if(Temp2!=nullptr)
                {
                    Temp2->Left = LeftMost->Right;
                }
                else{
                    Temp1->Right = LeftMost->Right;
                }
                Temp1->Data = LeftMost->Data;
            }
        }

        Node<Type>* getRootTree()
        {
            return this->Root;
        }

        bool hasItem(const Type FindingItem)
        {
            Node<Type>* Temp = Root;
            while(Temp != nullptr)
            {
                if(Temp->Data == FindingItem)
                {
                    return true;
                }

                if(FindingItem < Temp->Data)
                {
                    Temp = Temp->Left;
                }
                else{
                    Temp = Temp->Right;
                }
            }
            return false;
        }

        void clear(Node<Type>* RootTree)
        {
            if(RootTree == nullptr)
            {
                return;
            }

            clear(RootTree->Left);
            clear(RootTree->Right);
            delete RootTree;
            RootTree = nullptr;
        }

        void print(Node<Type>* RootTree, UINT32 Space=0, UINT16 t=0)
        {
            UINT32 Count = 3;
            if(RootTree == nullptr)
                return;
            Space += Count;

            print(RootTree->Right, Space, 1);

            for(UINT32 i=Count; i<Space; i++)
            {
                std::cout <<" ";
            }

            if(t==1) //Right node
            {
                std::cout <<"/"<< RootTree->Data <<"\n";
            }
            else if(t==2) //Left node
            {
                std::cout <<"\\"<< RootTree->Data <<"\n";
            }
            else{ //Root node
                std::cout << RootTree->Data <<"\n";
            }
            print(RootTree->Left, Space, 2);
        }
    };
//-------------------------------------------------------------------------------------------------------------------------------
    template<typename Type>
    class Matrix
    {
    private:
        UINT_T col;     // i - column
        UINT_T row;     // j - row
        Type** Array;
    public:
        Matrix()
        {
            col = 1;
            row = 1;
            Array = new Type*[col];
            Array[col] = new Type[row];
            Array[col-1][row-1] = 0;
        }
        Matrix(UINT_T col, UINT_T row)
        {
            this->col = col;
            this->row = row;

            Array = new Type*[col];
            for(UINT_T i=0; i<col; i++)
            {
                Array[i] = new Type[row];
            }
            for(UINT_T i=0; i<col; i++)
            {
                for(UINT_T j=0; j<row; j++)
                {
                    Array[i][j] = 0;
                }
            }

        }
        Matrix(const Matrix& OtherMatrix)
        {
            this->col = OtherMatrix.col;
            this->row = OtherMatrix.row;
            //this->Array = OtherMatrix.Array;

            Array = new Type*[this->col];
            for(UINT_T i=0; i<col; i++)
            {
                Array[i] = new Type[this->row];
            }

            for(UINT_T i=0; i<col; i++)
            {
                for(UINT_T j=0; j<row; j++)
                {
                    Array[i][j] = OtherMatrix.Array[i][j];
                }
            }
        }
        ~Matrix()
        {
            for(UINT_T i=0; i<col; i++)
            {
                if(Array[i] != nullptr)
                    delete[] Array[i];
            }
            if(Array != nullptr)
                delete[] Array;
        }

        void ScanArray()
        {
            for(UINT_T i=0; i<col; i++)
            {
                for(UINT_T j=0; j<row; j++)
                {
                    std::cout << "a["<<i+1<<"]["<<j+1<<"] = ";
                    std::cin >> Array[i][j];
                }
            }
        }

        void PrintArray()
        {
            for(UINT_T i=0; i<col; i++)
            {
                std::cout <<"( ";
                for(UINT_T j=0; j<row; j++)
                {
                    std::cout << Array[i][j] <<" ";
                }
                std::cout <<" )\n";
            }
        }
    };
//-------------------------------------------------------------------------------------------------------------------------------
}//END of Namespace SGM
#endif // SGM_HEADER_H_INCLUDED