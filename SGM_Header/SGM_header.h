/**
 *  Useful functions
 *  In this header used some standart functions of C++17
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
 *  14) GeneratePassword()
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
 *  13) Graph 
 *
 */
#ifndef SGM_HEADER_H_INCLUDED
#define SGM_HEADER_H_INCLUDED
#define DEBUG 1

// C# or Java style class modifiers
#define Public public:
#define Private private:
#define Protected protected:

// C# style
#define var auto
#define foreach for

#include <iostream>
#include <iomanip>
#include <string>
#include <cmath>
#include <list>
#include <stack>
#include <vector>
#include <stdlib.h>
//#include <utility>

namespace SGM 
{	
	//'using' syntax is new standart of C++11
    using uint64 = unsigned long long int;
	using uint32 = unsigned long int;
	using ushort = unsigned short int;
	using int64 = long long int;
	using int32 = long int;
	using uchar = unsigned char;
	using ldouble = long double;

	using std::string;
	using std::cout;
	using std::cin;
	using std::cerr;
	using std::swap;
	using std::runtime_error;
	using std::list;
	using std::vector;
	using std::stack;

//-------------------------------------------------------------------------------------------------------------------------------
//----------------------------------------------  FUNCTIONS  --------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------------------------
    /** Factorial  n! */
    uint64 Factorial(size_t inputNum)
    {
        uint64 result=1;
        for(uint32 i=1; i<=inputNum; i++)
        {
            result*=i;
        }
        return result;
    }

    /** Count digits in the integer number */
    template<typename int_t>
    uint64 CountDigits(int_t inputNum)
    {
        uint32 Count=0;
        while(inputNum !=0)
        {
            inputNum /= 10;
            Count++;
        }
        return Count;
    }


    /** Divide inputNum by digits , Its will work only for integer number!!! */
    template<typename int_t>
    void DivideNumsByDig(uint64 inputNumber, int_t* array)
    {
        uint32 len = CountDigits(inputNumber);

        for(uint32 i=0; i<len; i++)
        {
            if(i==0)
            {
                array[i] = inputNumber/pow(10,len-1);
                continue;
            }
            array[i] = inputNumber/pow(10,(len-1)-i);
            array[i] %= 10;
        }
    }

    /** Divide inputNum(type string) by digits numbers and convert to integer type */
    template<typename int_t>
    void DivideNumsByDig(string inputNumber, int_t* array)
    {
        uint32 len = inputNumber.length();
        for(uint32 i=0; i<len; i++)
        {
            char temp = inputNumber[i];
            array[i] = atoi(&temp); // char to int
        }
    }

    /** Get maximum of number in the array */
    template <typename Type>
    Type GetMax(const Type* array, uint32 size)
    {
        Type Max = array[0];

        for(uint32 i=0; i<size; i++)
        {
            if(Max<array[i])
                Max = array[i];
        }
        return Max;
    }

    /** Get maximum of number in the range: StartIterator...EndIterator */
    template <typename FowardIterator>
    FowardIterator GetMax(FowardIterator* start, FowardIterator* end)
    {
        FowardIterator* Max = start;
        uint32 Size = end - start;

        for(uint32 i=0; i<=Size; i++)
        {
            if(*(start+i) > *Max)
            *Max = *(start+i);
        }
        return *Max;
    }

    /** Get maximum of 2 numbers */
    template <typename Type>
    inline Type GetMax(Type firstNumber, Type secondNumber)
    {
        return (firstNumber>=secondNumber ? firstNumber:secondNumber);
    }

    /** Get minimum of number in the array */
    template <typename Type>
    Type GetMin(const Type* array, uint32 size)
    {
        Type Min = array[0];

        for(uint32 i=0; i<size; i++)
        {
            if(Min>array[i])
                Min = array[i];
        }
        return Min;
    }

    /** Get minimum of number in the range: StartIterator...EndIterator */
    template <typename FowardIterator>
    FowardIterator GetMin(FowardIterator* start, FowardIterator* end)
    {
        FowardIterator* Min = start;
        uint32 Size = end - start;

        for(uint32 i=0; i<=Size; i++)
        {
            if(*(start+i) < *Min)
                *Min = *(start+i);
        }
        return *Min;
    }

    /** Get minimum of 2 numbers */
    template <typename Type>
    inline Type GetMin(Type firstNumber, Type secondNumber)
    {
        return (firstNumber<=secondNumber ? firstNumber:secondNumber);
    }


    /** Sort 1D array by increasing order(Algorithm Bubble Sort, complex O(N^2)) */
    template <typename ArrayType>
    void SortInc_B(ArrayType* array, uint32 size)
    {
        ArrayType Temp; // temporary variable

        //Algorithm "Buble sort"
        for(uint32 i=0; i<size; i++)
        {
            for(uint32 j=0; j<(size-i-1); j++)
            {
                if(array[j]>array[j+1])
                {
                    Temp = array[j];
                    array[j] = array[j+1];
                    array[j+1]= Temp;
                }
            }
        }
    }

    /** Sort 1D array by decreasing order(Algorithm Bubble Sort, complex O(N^2)) */
    template <typename ArrayType>
    void SortDec_B(ArrayType* array, uint32 size) // Sort array by decreasing order
    {
        ArrayType Temp; // temporary variable

        //Algorithm "Buble sort"
        for(uint32 i=0; i<size; i++)
        {
            for(uint32 j=0; j<(size-i-1); j++)
            {
                if(array[j]<array[j+1])
                {
                    Temp = array[j];
                    array[j] = array[j+1];
                    array[j+1]= Temp;
                }
            }
        }
    }


    /** Input whatever array in the Console */
    template <typename Type>
    inline void ScanArray(Type* array, uint32 size)
    {
        for(uint32 i=0; i<size; i++)
        {
            cin >> array[i];
        }
    }

    /** Print whatever array to the display */
    template <typename Type>
    inline void PrintArray(const Type* array, uint32 size)
    {
        for(uint32 i=0; i<size; i++)
        {
            cout << array[i] << std::setw(4);
        }
    }

    /** Recursive Binary Search complex O(log(n)) */
    template <typename Type>
    int BinarySearch_Recursive(const Type* array, uint32 left, uint32 right, const Type findElement)
    {
        //int32 Left = 1;
        //int32 Right = Size;
        int Middle = (left + right)/2;

        if(right >= left)
        {           
            Middle = (left + right)/2;

            if(array[Middle] == findElement)
            {
                return Middle;                
            }
            if(array[Middle] > findElement)
            {
                return Recursive_BinarySearch(array, left, Middle - 1, findElement);
            }
            else //if(Array[Middle] < FindElement)
            {
                return Recursive_BinarySearch(array, Middle + 1, right, findElement);
            }
        }
		return -1;
    }

    /** Simple Binary Search complex O(log(n)), if element was found in the sorted array, then function returned position of element in sorted array, otherwise function will return -1 */
    template <typename Type>
    int BinarySearch(const Type* array, uint32 size, Type findElement)
    {
        int32 low = 1;
        int32 high = size;
        bool Founded = false;

        while(low<=high)
        {
            int32 mid = (low+high)/2;

            if(array[mid]>findElement)
            {
                high = mid-1;
            }
            else if(array[mid]<findElement)
            {
                low = mid+1;
            }
            else if(array[mid]==findElement)
            {
                return mid;                
            }
        }
        return -1;
    }

    /** Function calculate the greatest common divisor(GCD) two natural numbers */
    template<typename int_t>
    int_t GCD(int_t firstNumber, int_t secondNumber)
    {
        while(firstNumber*secondNumber>0)
        {
            if(firstNumber>=secondNumber)
            {
                firstNumber %= secondNumber;
            }
            else{
                secondNumber %= firstNumber;
            }
        }
        return firstNumber+secondNumber; //ÍÎÄ
    }

    /** Function calculate the least common multiples(LCM) two natural numbers */
    template<typename int_t>
    int_t LCM(int_t firstNumber, int_t secondNumber)
    {
        return (firstNumber*secondNumber) / GCD(firstNumber, secondNumber); //ÍÎÊ (íàèìåíüøåå îáùåå êðàòíîå)
    }

	/**
	* @param PassLength - Password Length
	* @param NumberOfUpperAlpha - Number of Upper Alphabets must be in the password
	* @param NumberOfLowerAlpha - Number of Lower Alphabets must be in the password
	* @param NumberOfDigit - Number of Digits must be in the password
	* @param HaveTwoIdenChar - Have two identical consecutive characters in the password, by default its activated
	* @return Function will return generated password which satisfying the requirements
	*/
	string GeneratePassword(size_t passLength, size_t numberOfUpperAlpha, size_t numberOfLowerAlpha, size_t numberOfDigit, bool haveTwoIdenChar = false)
	{
		srand((size_t)time(false));
		const string Chars[] = { "ABCDEFGHIJKLMNOPQARTUVWXYZ", "abcdefghijklmnopqrstuvwxyz", "0123456789" };
		const size_t UpperAlphLen = Chars[0].length();
		const size_t LowerAlphLen = Chars[1].length();
		const size_t NumbersLen = Chars[2].length();
		string Password;

		if ((numberOfUpperAlpha + numberOfLowerAlpha + numberOfDigit) > passLength)
		{
			throw runtime_error("Argument ERROR: Number of upper alphabets and Number of lower alphabets and Number of digits more than Password length \n");
		}

		for (size_t i = 1; i <= numberOfUpperAlpha; i++)
		{
			Password.push_back(Chars[0][0 + rand() % UpperAlphLen]);
		}
		for (size_t i = 1; i <= numberOfLowerAlpha; i++)
		{
			Password.push_back(Chars[1][0 + rand() % LowerAlphLen]);
		}
		for (size_t i = 1; i <= numberOfDigit; i++)
		{
			Password.push_back(Chars[2][0 + rand() % NumbersLen]);
		}

		for (size_t i = 1; i <= passLength; i++)
		{
			swap(Password[0 + rand() % passLength], Password[0 + rand() % passLength]);
		}

		if (!haveTwoIdenChar)
		{
			for (size_t i = 0; i<passLength; i++)
			{
				if (i != passLength - 1 && Password[i] == Password[i + 1])
				{
					string Temp = Chars[0 + rand() % 2];
					Password[i] = Temp[0 + rand() % Temp.length()];
				}
			}
		}

		return Password;
	}

//-------------------------------------------------------------------------------------------------------------------------------
//----------------------------------------------  CLASSES  ----------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------------------------
	struct HashEntry
	{
		string Value;

		HashEntry(string Value)
		{
			this->Value = Value;
		}
	};
//-------------------------------------------------------------------------------------------------------------------------------
    /** class Hash Table */
	class HashTable
    {
    private:
		const int TABLE_SIZE = 127; // Èñïîëüçîâàíî ïðîñòî ÷èñëî
		HashEntry **Table; //Òàáëèöà õåøîâ

		int64 HashFunction(string Value)
		{
			int64 Sum = 0; //ASCII sum
			for (var i : Value)
			{
				Sum += i;
			}
			return Sum % TABLE_SIZE;
		}
    public:
		HashTable()
		{
			Table = new HashEntry*[TABLE_SIZE];
			for (int i = 0; i < TABLE_SIZE; i++)
			{
				Table[i] = nullptr;
			}
		}
		HashTable(const HashTable& OtherHashTable)
		{
			//this->TABLE_SIZE = OtherHashTable.TABLE_SIZE;
			this->Table = OtherHashTable.Table;
			for (int i = 0; i < this->TABLE_SIZE; i++)
			{
				this->Table[i] = OtherHashTable.Table[i];
			}
		}
        ~HashTable()
        {
			for (int i = 0; i < TABLE_SIZE; i++)
			{
				delete Table[i];
			}
			delete Table;
        }

		void insert(string Value)
		{
			int64 HashNumber = HashFunction(Value);
			if (Table[HashNumber] == nullptr)
			{
				Table[HashNumber] = new HashEntry(Value);
			}
		}

		void remove(string Value)
		{
			int64 HashNumber = HashFunction(Value);
			if (Table[HashNumber] != nullptr)
			{
				Table[HashNumber] = nullptr;
			}
		}

		HashEntry* find(string Value)
		{
			int64 HashNumber = HashFunction(Value);
			if (Table[HashNumber] == nullptr)
			{
				return nullptr;
			}
			else{
				return Table[HashNumber];
			}
		}

		int64 getHash(string Value)
		{
			int64 HashNumber = HashFunction(Value);
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
            cout <<"X="<< X <<" Y="<< Y <<" Z="<< Z <<"\n";
        }
    };

//-------------------------------------------------------------------------------------------------------------------------------

    template <class T>
    class Determinant
    {
    private:
        T det; // element of determinant
        T result;
        uint32 size_det;

        void GetMinor(T** (&Array), T** p, size_t i, size_t j, size_t m)
        {
            size_t ki=0, kj=0, di=0, dj=0;

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

        void Display(T** (&Array), uint32 Size)
        {
            Size = size_det;
            Array = det;

            cout << "\n\n" << "|";
            for(uint32 i=0; i<size_det; i++)
            {
                uint32 j=0;
                for(; j<size_det; j++)
                {
                    cout << std::setw(4) << det[i][j] << std::setw(4);
                }
                cout << "|" << std::endl;

                if((i != size_det-1) && (j != size_det-1)){
                    cout << "|";
                }
            }
            cout << "\n\n";
        }

        T Compute(T** (&Array), uint32 Size) // Recruisive computing determinant
        {
            Size = size_det;
            Array = det;

            T **p = nullptr;

            int k = 1; //(-1)^(i+j)
            uint32 i=0,j=0;

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
                    result = result + (k*(det[i][0])*ComputeDeterminant(p,(size_det-1)));
                    k = -k;
                }
            }
        }

        void ShowResult(T** (&array), uint32 Size)
        {
            cout << "D=" << result;
        }
    };
//-------------------------------------------------------------------------------------------------------------------------------
    /** class Big Integer */
    class BigInteger final
    {
    private:
        size_t* ArrayDigits;
        size_t SIZE;
        bool isCorrectNumber = true;

    public:
        const static size_t ONE = 1;
        const static size_t ZERO = 0;

        BigInteger()
        {
            SIZE = 1;
            ArrayDigits = new size_t[SIZE];
            ArrayDigits[0] = this->ZERO;

        }

        BigInteger(const BigInteger *A)
        {
            SIZE = A->SIZE;
            ArrayDigits = new size_t[SIZE];
            for(uint64 i=0; i<SIZE; i++)
            {
                ArrayDigits[i] = A->ArrayDigits[i];
            }
        }

        BigInteger(string BigIntegerNumber)
        {
            SIZE = BigIntegerNumber.length();
            ArrayDigits = new size_t[SIZE];

            for(int i=SIZE-1; i>=0; i--)
            {
                char temp = BigIntegerNumber[(SIZE-1)-i];

                try
                {
                    if(!isdigit(temp))
                    {
                        isCorrectNumber = false;
                        throw runtime_error("\nConvert Error Big Integer \n");
                    }
                    else{
                        size_t num = atoi(&temp); // char to int

                        if(num>=10)
                        {
                            num /= 10;
                        }
                        ArrayDigits[i] = num;
                    }
                }
                catch(runtime_error& e)
                {
                    cerr << e.what();
                    break;
                }
            }
        }

        ~BigInteger()
        {
            if(ArrayDigits != nullptr)
            {
                cout<<"Destructor called\n";
                delete[] ArrayDigits;
            }
        }

        string ToString()
        {
            string str;
            string temp;

            for(int64 i=SIZE-1; i>=0 && isCorrectNumber; i--)
            {
                temp = std::to_string(ArrayDigits[i]);
                str.push_back(temp[0]);
                cout << str << "\n";
            }
            return str;
        }

        uint64 getSize()
        {
            return SIZE;
        }

        BigInteger& operator=(const BigInteger*)
        {
            return *this;
        }

        BigInteger& operator=(string BigIntegerNumber)
        {
            SIZE = BigIntegerNumber.length();
            ArrayDigits = new size_t[SIZE];

            for(int i=SIZE-1; i>=0; i--)
            {
                char temp = BigIntegerNumber[(SIZE-1)-i];

                try
                {
                    if(!isdigit(temp))
                    {
                        isCorrectNumber = false;
                        throw runtime_error("\nConvert Error Big Integer \n");
                    }
                    else{
                        size_t num = atoi(&temp); // char to int

                        if(num>=10)
                        {
                            num /= 10;
                        }
                        ArrayDigits[i] = num;
                    }
                }
                catch(runtime_error& e)
                {
                    cerr << e.what();
                    break;
                }
            }
            return *this;
        }

        BigInteger& operator+(const BigInteger &A)
        {
            BigInteger* res = new BigInteger(); //res - result
            size_t carry = 0;

            if(SIZE == A.SIZE)
            {
                for(uint64 i=0; i<SIZE; i++)
                {
                    size_t sum = ArrayDigits[i] + A.ArrayDigits[i] + carry;
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
            for(uint64 i=0; i<SIZE; i++)
            {
                cout<<ArrayDigits[i]<<"\n";
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
        uint64 Size; //Размер стека
        uint64 Count; //Количество добавленных элементов в стеке
        Type* Array; //Главный массив
        Type* NewArray = nullptr; //Вспомогательный массив

    public:
        Stack(uint64 MaxSize = 0) 
        {
            Size = MaxSize;
            Count = 0;
            if(MaxSize==0)
                Array = new Type[++MaxSize]{0};
            else
                Array = new Type[MaxSize]{0};
        }

        Stack(const Stack<Type> &OtherStack)
        {
            Size = OtherStack.Size;
            Count = OtherStack.Count;
            Array = new Type[Size];

            for(uint64 i=0; i<Size; i++)
            {
                Array[i] = OtherStack.Array[i];
            }
        }

        ~Stack()
        {
            if(Array != nullptr)
                delete[] Array;
            if(NewArray != nullptr)
                delete[] NewArray;
        }

        void push_back(const Type Value)
        {
            NewArray = new Type[++Size]{0};
            for(uint64 i=0; i<Size; i++)
            {
                NewArray[i] = Array[i];
            }
            delete[] Array;
            Array = NewArray;

            Array[Count++] = Value;
        }

        void pop_back()
        {
            try
            {
                if(Size==0 || Count==0)
                    throw runtime_error("\nERROR: the Size of Stack is NULL \n");

                NewArray = new Type[Size]{0};
                for(uint64 i=0; i<Size; i++)
                {
                    NewArray[i] = Array[i];
                }
                delete[] Array;
                Array = NewArray;

                Size--;
                Count--;
                Array[Count] = 0;
            }
            catch(runtime_error& e)
            {
                cerr << e.what();
            }
        }

        Type back() 
        {
            if(Count!=0)
                return Array[Count-1];
            else
                return 0;
        }

        uint64 size() const 
        {
            return Count;
        }
       
        bool empty() 
        {
            if(Size>=1 && Count>0)
                return false;
            else
                return true;
        }
        void clear() 
        {
            for(uint64 i=0; i<Size; i++)
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

            for(uint64 i=0; i<Size; i++)
            {
                this->Array[i] = OtherStack.Array[i];
            }
            return *this;
        }

        #if DEBUG >= 1
        void print()
        {
            for(uint64 i=0; i<Size; i++)
            {
                cout << Array[i] <<" ";
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
        uint64 Size;    //Размер Очереди
        uint64 Begin;   //Начало Очереди
        uint64 End;     //Конец Очереди
        Type* Array;    //Главный массив
        Type* NewArray = nullptr; //Вспомогательный массив

    public:
        Queue(uint64 MaxSize = 0) 
        {
            Size = MaxSize;
            Begin = 0;
            End = 0;
            if(MaxSize==0)
                Array = new Type[++MaxSize]{0};
            else
                Array = new Type[MaxSize]{0};
        }

        Queue(const Queue<Type> &OtherQueue)
        {
            Size = OtherQueue.Size;
            Begin = OtherQueue.Begin;
            End = OtherQueue.End;
            Array = new Type[Size];

            for(uint64 i=0; i<Size; i++)
            {
                Array[i] = OtherQueue.Array[i];
            }
        }

        ~Queue() 
        {
            if(Array != nullptr)
                delete[] Array;
            if(NewArray != nullptr)
                delete[] NewArray;
        }

        void push_back(const Type Value) 
        {
            NewArray = new Type[++Size]{0};
            for(uint64 i=0; i<Size; i++)
            {
                NewArray[i] = Array[i];
            }
            delete[] Array;
            Array = NewArray;

            End = Size - 1;
            Array[End] = Value;
        }

        void pop_front() 
        {
            try
            {
                if(Size==0 || End<0 || Begin>End)
                    throw runtime_error("\nERROR: the Size of Queue is NULL \n");

                NewArray = new Type[Size]{0};
                for(uint64 i=0; i<Size; i++)
                {
                    if(i==Size-1)
                    {
                        NewArray[i] = 0; 
                    }
                    else{
                        NewArray[i] = Array[i+1];
                    }
                }

                delete[] Array;
                Array = NewArray;

                Size--;
                End = Size - 1;
                //Begin = 0;
            }
            catch(runtime_error& e)
            {
                cerr << e.what();
            }

        }

        Type front() 
        {
            if(End>=Begin)
                return Array[Begin];
            else
                return 0;
        }

        Type back()
        {
            if(End>=Begin)
                return Array[End];
            else
                return 0;
        }

        uint64 size() const 
        {
            return Size;
        }

        bool empty() 
        {
            if(Size>=1 && End>=Begin)
                return false;
            else
                return true;
        }

        void clear() 
        {
            for(uint64 i=0; i<Size; i++)
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

            for(uint64 i=0; i<Size; i++)
            {
                this->Array[i] = OtherQueue.Array[i];
            }
            return *this;
        }

        #if DEBUG >= 1
        void print()
        {
            cout << "Begin: "<<Begin<<" End: "<<End<<"\n";
            for(uint64 i=0; i<Size; i++)
            {
                cout << Array[i] <<" ";
            }
            //cout <<"\n"<< Array[End+1] <<" \n";
        }
        #endif // DEBUG
    };

//-------------------------------------------------------------------------------------------------------------------------------
    /** Template of Class Deque, Deque is structure data that realized the interface of "Double-Ended-Queue"  */
    template<class Type>
    class Deque
    {
    private:
        uint64 Size;                //Размер Дека
        uint64 Begin;               //Начало Дека (Iterator begin)
        uint64 End;                 //Конец Дека (Iterator end)
        Type* Array;                //Главный массив
        Type* NewArray = nullptr;   //Вспомогательный массив

        typedef Type* iterator;
        typedef Type& reference;
        typedef const Type* const_iterator;
        typedef const Type& const_reference;

    public:
        Deque(uint64 MaxSize = 0) 
        {
            Size = MaxSize;
            Begin = 0;
            End = 0;
            
            if(MaxSize==0)
                Array = new Type[++MaxSize]{0};
            else
                Array = new Type[MaxSize]{0};
        }

        Deque(const Deque<Type> &OtherDeque) 
        {
            Size = OtherDeque.Size;
            Begin = OtherDeque.Begin;
            End = OtherDeque.End;
            Array = new Type[Size];

            for(uint64 i=0; i<Size; i++)
            {
                Array[i] = OtherDeque.Array[i];
            }
        }

        ~Deque() 
        {
            if(Array != nullptr)
                delete[] Array;
            if(NewArray != nullptr)
                delete[] NewArray;
        }

        void push_back(const Type Value) 
        {
            NewArray = new Type[++Size]{0};
            for(uint64 i=0; i<Size; i++)
            {
                NewArray[i] = Array[i];
            }
            delete[] Array;
            Array = NewArray;

            End = Size - 1;
            Array[End] = Value;
        }

        void push_front(const Type Value) 
        {
            NewArray = new Type[++Size]{0};
            for(uint64 i=0; i<Size; i++)
            {
                if(i!=Size-1)
                {
                    NewArray[i+1] = Array[i]; 
                }
            }
            delete[] Array;
            Array = NewArray;

            End = Size - 1;
            Array[Begin] = Value;                                     
        }

        void pop_front() 
        {
            try
            {
                if(Size==0 || End<0 || Begin>End)
                    throw runtime_error("\nERROR: the Size of Queue is NULL \n");

                NewArray = new Type[Size]{0};
                for(uint64 i=0; i<Size; i++)
                {
                    if(i==Size-1)
                    {
                        NewArray[i] = 0; //Îáíóëèì êîíå÷íûé ýëåìåíò â ìàññèâå
                    }
                    else{
                        NewArray[i] = Array[i+1]; //Ñäâèã ýëåìåíòîâ ìàññèâà â âëåâî(<--)
                    }
                }

                delete[] Array;
                Array = NewArray;

                Size--;
                End = Size - 1;
                //Begin = 0;
            }
            catch(runtime_error& e)
            {
                cerr << e.what();
            }
        }

        void assign(const Type Value, const size_t Pos)
        {
            try
            {
                if(Pos>=Size)
                {
                    throw runtime_error("\nERROR: the index of element is largest than Deque size \n");
                }
                this->Array[Pos] = Value;
            }
            catch(runtime_error& e)
            {
                cerr << e.what();
                return;
            }
        }

        void pop_back() 
        {
            try
            {
                if(Size==0 || End<0 || Begin>End)
                    throw runtime_error("\nERROR: the Size of Deque is NULL \n");

                NewArray = new Type[Size]{0};
                for(uint64 i=0; i<Size; i++)
                {
                    NewArray[i] = Array[i];
                }
                delete[] Array;
                Array = NewArray;

                Size--;
                End = Size - 1;              
            }
            catch(runtime_error& e)
            {
                cerr << e.what();
            }
        }

        const_reference front() const
        {
            //if(End>=Begin)
                return *(this->Array+this->Begin);
        }

        const_reference back() const 
        {
            //if(End>=Begin)
                return *(this->Array+this->End);
        }

        uint64 size() const 
        {
            return Size;
        }

        bool empty() 
        {
            if(Size>=1 && End>=Begin)
                return false;
            else
                return true;
        }

        void resize(size_t NewSize)
        {
            try
            {
                if(NewSize<=0)
                    throw runtime_error("\nERROR: the Size of Deque is least than ZERO \n");

                NewArray = new Type[NewSize]{0};
                if(NewSize>=Size)
                {
                    for(uint64 i=0; i<Size; i++)
                    {
                        NewArray[i] = Array[i];
                    }
                }
                else{
                    for(uint64 i=0; i<NewSize; i++)
                    {
                        NewArray[i] = Array[i];
                    }
                }
                delete[] Array;
                Array = NewArray;

                Size = NewSize;
                End = Size - 1;

            }
            catch(runtime_error& e)
            {
                cerr << e.what();
            }
        }

        void clear() 
        {
            for(uint64 i=0; i<Size; i++)
            {
                Array[i] = 0;
            }
            Begin = 0;
            End = 0;
            Size = 0;            
        }

        const_iterator begin() const
        {
            return this->Array + this->Begin;
        }
        const_iterator end() const 
        {
            return this->Array + this->End;
        }

        const_reference at(size_t ElementIndex) const
        {
            try
            {             
                if(ElementIndex>=Size)
                {
                    throw runtime_error("\nERROR: the index of element is largest than Deque size \n");
                }
                return *(this->Array+ElementIndex);
            }
            catch(runtime_error& e)
            {
                cerr << e.what();
                return 0;
            }
        }

        const_reference operator[](size_t ElementIndex) const
        {           
            return *(this->Array+ElementIndex);
        }

        Deque<Type>& operator=(const Deque<Type>& OtherDeque)
        {
            this->Size = OtherDeque.Size;
            this->Begin = OtherDeque.Begin;
            this->End = OtherDeque.End;

            Array = new Type[this->Size];
            for(uint64 i=0; i<this->Size; i++)
            {
                Array[i] = OtherDeque.Array[i];
            }

            return *this;
        }     


        #if DEBUG >= 1
        void print()
        {
            cout << "Begin: "<<Begin<<" End: "<<End<<"\n";
            for(uint64 i=0; i<Size; i++)
            {
                cout << Array[i] <<" ";
            }
            //cout <<"\n"<< Array[End+1] <<" \n";
        }
        #endif // DEBUG
    };
//-------------------------------------------------------------------------------------------------------------------------------
    /** Template struct of Node, Node is structure data which saved itself data and address of next and previous data */
    template<class Type>
    struct BinaryNode
    {
        Type Data;          //data
        BinaryNode<Type>* Right;  //pointer to right node (used for binary tree)
        BinaryNode<Type>* Left;   //pointer to left node (used for binary tree)

        BinaryNode(): Right(nullptr), Left(nullptr) {} //Default constructor
        BinaryNode(const Type& setData): Data(setData),Right(nullptr), Left(nullptr) {}
        BinaryNode(const BinaryNode<Type>& OtherNode) //Copying constructor
        {
            this->Data = OtherNode.Data;
            this->Left = OtherNode.Left;
            this->Right = OtherNode.Right;
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
        BinaryNode<Type>* Head;     //Указатель на начало списка
        typedef BinaryNode<Type>* iterator;
        typedef const BinaryNode<Type>* const_iterator;
        typedef Type& reference;
        typedef const Type& const_reference;

    public:
        ForwardList(): Head(nullptr) {}
        ForwardList(const ForwardList<Type>& OtherForwardList)
        {
            this->Head = OtherForwardList.Head;
        }
        ~ForwardList()
        {
            while(Head!=nullptr)  //Пока адрес не пусто
            {
                BinaryNode<Type>* Temp = Head->Right;  //Временная переменная для хранения адреса следующего элемента
                delete Head;                     //Освобождаем адрес обозначающий начало
                Head = Temp;                     //Меняем адрес на следующий
            }
        }

		//Метод добавления элементов в список
        void push_front(const Type Value) 
        {
            BinaryNode<Type>* Temp = new BinaryNode<Type>(); //При каждом вызове выделяется память
            Temp->Data = Value;  //Записываем Value в элемент структуры  Node
            Temp->Right = Head;  //Указываем, что след. элемент это объект по адресу Head
            Head = Temp;         //Указываем, что последний активный элемент это только что введенный
        }

		//Метод добавления элементов в список
        void push_back(const Type Value) 
        {
            BinaryNode<Type>* Temp = new BinaryNode<Type>();
            BinaryNode<Type>* NewNode = new BinaryNode<Type>(Value);
            Temp = Head;
            while(Temp->Right != nullptr)
            {
                Temp = Temp->Right;
            }
            Temp->Right = NewNode;
        }

        void assign(const Type Value, const size_t Pos)
        {
            if(Head!=nullptr && Pos <= this->size()) //Делаем проверку на то что список не пуст и Pos не превышает число его элементов
            {
                BinaryNode<Type>* Temp = Head;
                for(uint64 i=1; i<=Pos; i++)
                {
                    Temp = Temp->Right; 
                }
                Temp->Data = Value;
            }
        }

        void insert(const Type Value, const size_t Pos)
        {
            if(Head!=nullptr && Pos!=1 && Pos < this->size() ) 
            {
                BinaryNode<Type>* BeforeNode = Head;
                BinaryNode<Type>* AfterNode = Head;
                BinaryNode<Type>* NewNode = new BinaryNode<Type>(Value);

                for(uint64 i=1; i<Pos-1; i++)
                {
                    BeforeNode = BeforeNode->Right;
                }
                for(uint64 i=1; i<Pos; i++)
                {
                    AfterNode = AfterNode->Right;
                }
                BeforeNode->Right = NewNode;
                NewNode->Right = AfterNode;             
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

		//Метод удаления элемента из списке
        void erase(const Type Element) 
        {
            BinaryNode<Type>* ElementNode = Head;
            BinaryNode<Type>* Temp = Head;

            
            while(ElementNode->Right!=nullptr && ElementNode->Data != Element)
            {
                ElementNode = ElementNode->Right;
            }

            
            if(Temp==ElementNode)
            {
                Temp->Right = ElementNode->Right;
            }
            else{
                while(Temp->Right != ElementNode)
                {
                    Temp = Temp->Right;
                }
                Temp->Right = ElementNode->Right;
            }
        }

		//Удаление одного звенья из начало списка
        void pop_front() 
        {
            Head = Head->Right;
        }

		//Удаление N звен из начало списка
        void pop_front_nodes(const size_t N) 
        {
            if(Head!=nullptr  &&  N <= this->size()) 
            {
                for(uint64 i=1; i<=N; i++)
                {
                    this->pop_front();
                }
            }
        }

		//Удаление одного звенья из конца списка
        void pop_back() 
        {
            BinaryNode<Type>* Temp = Head;
            BinaryNode<Type>* LastElement = nullptr;

            while(Temp->Right!=nullptr) 
            {
                Temp = Temp->Right;
            }
            LastElement = Temp;
            LastElement->Data = Temp->Data; 
            this->erase(LastElement->Data); 
        }

		//Удаление N звен из конца списка
        void pop_back_nodes(const size_t N) 
        {
            if(Head!=nullptr  &&  N <= this->size()) 
            {
                for(uint64 i=1; i<=N; i++)
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
            BinaryNode<Type>* Temp = Head;
            while(Temp->Right!=nullptr)
            {
                Temp = Temp->Right;
            }
            return Temp->Data;
        }

        const_iterator begin() const
        {
            return Head;
        }

        const_iterator end() const
        {
            BinaryNode<Type>* Temp = Head;
            while(Temp->Right!=nullptr)
            {
                Temp = Temp->Right;
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
            while (Head!=nullptr) 
            {
                BinaryNode<Type>* Temp = Head->Right;  
                delete Head;                    
                Head = Temp;                    
            }
        }

        bool empty()
        {
            if(this->size())
                return false;
            else
                return true;
        }

        uint64 size()
        {
            uint64 Count = 0;
            BinaryNode<Type>* Temp = Head;
            while(Temp!=nullptr)
            {
                Count++;
                Temp = Temp->Right;
            }
            return Count;
        }

        #if DEBUG >= 1
        void print()
        {
            BinaryNode<Type>* Temp = Head; 
            while(Temp!=nullptr) 
            {
                cout << Temp->Data <<" -> ";
                Temp = Temp->Right; 
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
        BinaryNode<Type>* Root;

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
            BinaryNode<Type>* Temp1 = new BinaryNode<Type>();
            if(Root==nullptr)
            {
                Root = Temp1;
                Root->Data = Item;
                return;
            }

            Temp1 = Root;
            BinaryNode<Type>* Temp2 = nullptr;

            while(Temp1 != nullptr) 
            {
                if(Temp1->Data == Item)
                {                  
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

            BinaryNode<Type>* NewNode = new BinaryNode<Type>(Item);
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
            BinaryNode<Type>* Temp1 = Root;
            BinaryNode<Type>* Temp2 = nullptr;

            while(Temp1 != nullptr)
            {
                if(Temp1->Data == Item)
                {                  
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
                BinaryNode<Type>* LeftMost = Temp1->Right;
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

        BinaryNode<Type>* getRootTree()
        {
            return this->Root;
        }

        bool hasItem(const Type FindingItem)
        {
            BinaryNode<Type>* Temp = Root;
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

        void clear(BinaryNode<Type>* RootTree)
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

        void print(BinaryNode<Type>* RootTree, uint32 Space=0, size_t t=0)
        {
            uint32 Count = 3;
            if(RootTree == nullptr)
                return;
            Space += Count;

            print(RootTree->Right, Space, 1);

            for(uint32 i=Count; i<Space; i++)
            {
                cout <<" ";
            }

            if(t==1) //Right node
            {
                cout <<"/"<< RootTree->Data <<"\n";
            }
            else if(t==2) //Left node
            {
                cout <<"\\"<< RootTree->Data <<"\n";
            }
            else{ //Root node
                cout << RootTree->Data <<"\n";
            }
            print(RootTree->Left, Space, 2);
        }
    };
//-------------------------------------------------------------------------------------------------------------------------------
    template<typename Type>
    class Matrix
    {
    private:
        size_t col;     // i - column
        size_t row;     // j - row
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
        Matrix(size_t col, size_t row)
        {
            this->col = col;
            this->row = row;

            Array = new Type*[col];
            for(size_t i=0; i<col; i++)
            {
                Array[i] = new Type[row];
            }
            for(size_t i=0; i<col; i++)
            {
                for(size_t j=0; j<row; j++)
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
            for(size_t i=0; i<col; i++)
            {
                Array[i] = new Type[this->row];
            }

            for(size_t i=0; i<col; i++)
            {
                for(size_t j=0; j<row; j++)
                {
                    Array[i][j] = OtherMatrix.Array[i][j];
                }
            }
        }
        ~Matrix()
        {
            for(size_t i=0; i<col; i++)
            {
                if(Array[i] != nullptr)
                    delete[] Array[i];
            }
            if(Array != nullptr)
                delete[] Array;
        }

        void Scan()
        {
            for(size_t i=0; i<col; i++)
            {
                for(size_t j=0; j<row; j++)
                {
                    cout << "a["<<i+1<<"]["<<j+1<<"] = ";
                    cin >> Array[i][j];
                }
            }
        }

        void Print()
        {
            for(size_t i=0; i<col; i++)
            {
                cout <<"( ";
                for(size_t j=0; j<row; j++)
                {
                    cout << Array[i][j] <<" ";
                }
                cout <<" )\n";
            }
        }
    };
//-------------------------------------------------------------------------------------------------------------------------------
	class Graph
	{
		Private int vertexes; //Number of vertexes
		Private list<int> *adjacencyList; //adjacency List
		
		Public Graph(int vertexes)
		{
			this->vertexes = vertexes;
			adjacencyList = new list<int>[vertexes];
		}

		Public Graph(const Graph& OtherGraph)
		{
			this->adjacencyList = OtherGraph.adjacencyList;
			this->vertexes = OtherGraph.vertexes;
		}

		Public void AddEdge(int vertex, int weight)
		{
			adjacencyList[vertex].push_back(weight);
		}

		Public void DepthFirstSearch(int sourceVertex)
		{
			//Initially mark all verices as not visited
			vector<bool> visited(this->vertexes, false);

			//Create a stack for DFS
			stack<int> Stack;

			//Push the current source node
			Stack.push(sourceVertex);

			while (!Stack.empty())
			{
				//Pop a vertex from stack and print it
				sourceVertex = Stack.top();
				Stack.pop();

				/* Stack may contain same vertex twice. So we need to print the popped item only
				   if it is not visited
				*/
				if (!visited[sourceVertex])
				{
					cout << sourceVertex << " -> ";
					visited[sourceVertex] = true;
				}

				/* Get all adjacent vertices of the popped vertex sourceVertex.
				   If a adjacent has not been visited, then push it to the stack
				*/
				foreach(var i : adjacencyList[sourceVertex])
				{
					if (!visited[i])
					{
						Stack.push(i);
					}
				}

			}
		}
	};
//-------------------------------------------------------------------------------------------------------------------------------
}//END of Namespace SGM
#endif // SGM_HEADER_H_INCLUDED
