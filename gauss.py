
def Gauss(matrix, column):
    '''
    Решение СЛАУ методом Гаусса. Принимаются два аргумента:\n
    matrix - строка строк, верхний уровень - число уравнений, нижний - набор коэффициентов одного уравенения;\n
    column - строка свободных членов
    '''
    start = [0 for i in range(len(column))]
    
    if len(matrix) != len(column):
        print("Number of the unkown don't equal number of the equation")
        return start
    
    # Число неизвестных/размер матрицы
    N = len(column)
    # Последовательность неизвестных
    # x = [i for i in range(N)]

    for i in range(N):
        # max_el = 0
        # max_id = None
        
        # for j in range(i, N):
        #     if abs(matrix[j][i]) >= max_el:
        #         max_el = matrix[j][i]
        #         max_id = j
        
        # matrix = matrix[:i] + matrix[max_id:max_id+1] + matrix[i:max_id] + matrix[max_id+1:]
        # column = column[:i] + column[max_id:max_id+1] + column[i:max_id] + column[max_id+1:]
        # x = x[:i] + x[max_id:max_id+1] + x[i:max_id] + x[max_id+1:]

        for j in range(i, N):
            cof = matrix[i][i]
            for k in range(i, N):
                matrix[j][k] /= cof# max_el
            column[j] /= cof# max_el

        for j in range(i + 1, N):
            cof = matrix[j][i]
            for k in range(i, N):
                matrix[j][k] -= cof * matrix[i][k]
            column[j] -= cof * column[i]
        
        # print(matrix)
        # print(column)
        # print(x)
    
    # print(matrix)
    # print(column)
    # print(x)

    res = []
    for i in range(N - 1, -1, -1):
        val = column[i]
        for j in range(N - 1, i, -1):
            val -= matrix[i][j] * res[N - 1 - j]
        res.append(val)
    res.reverse()

    return res


if __name__ == '__main__':
    print('Gauss test')
    matrix = [[1, 2, 3], [4, 5, 6], [18, 19, 10]]
    column = [4, 5, 6]
    print('Matrix is:')
    for i in range(len(matrix)):
        print(matrix[i])
    print('Is equal:')
    print(column)
    res = Gauss(matrix, column)
    print('Answer is: ')
    print(res)