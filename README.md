Найти точку максимума функции  
![eq](equations/func1.png)  
методом покоординатного спуска. Для одномерной минимизации использовать метод квадратичной интерполяции. Для поиска интервала унимодальности использовать алгоритм скользящего окна.

## Теоретическая часть
Задача отыскания максимума функции сводится к задаче нахождения минимума функции, противоположной по знаку заданной. Таким образом,
достаточно найти точку минимума функции  
![eq](equations/func2.png)  
С помощью метода покоординатного спуска функция минимизируется в направлениях, совпадающих с координатными осями и задаваемых
формулой (циклический покоординатный спуск):  
  
![eq](equations/coord1.png) где ![eq](equations/coord2.png) (единица стоит на j-ом месте), ![eq](equations/coord3.png).  
  
Минимизация осуществляется с помощью метода квадратичной интерполяции, который задает три пробные точки ![eq](equations/interp1.png), ![eq](equations/interp2.png) и ![eq](equations/interp3.png). Для нахождения точки ![eq](equations/interp2.png) задается шаг ![eq](equations/interp5.png) в положительном направлении от точки ![eq](equations/interp4.png), т.е. ![eq](equations/interp6.png) и если ![eq](equations/interp7.png), то  ![eq](equations/interp8.png), иначе  ![eq](equations/interp9.png)  
Строится квадратичный интерполяционый многочлен и находится его точка минимума по формуле:  
  
![eq](equations/interp10.png)  
Также находится точка ![eq](equations/interp11.png)  
Критерием окончания процесса является выполнение условий:  
![eq](equations/interp12.png),  ![eq](equations/interp13.png)  
Если условия окончания не выполняются и ![eq](equations/interp14.png), точка ![eq](equations/interp4.png) заменяется на точку ![eq](equations/interp15.png), в противном случае точка ![eq](equations/interp4.png) заменяется на ![eq](equations/interp16.png).  
  
Производить минимизацию возможно только на интервале унимодальности функции, который находится с помощью алгоритма скользящего окна: для выбранной исходной точки ![eq](equations/win1.png) и выбраного окна шириной ![eq](equations/win2.png) около точки ![eq](equations/win1.png) проверяется условие унимодальности:  
![eq](equations/win3.png)  
Если условие выполнено, то считается, что интервал унимодальности найден, в противном случае проверяется условие:  
![eq](equations/win4.png)  
Если оно выполнено, то ![eq](equations/win5.png), иначе ![eq](equations/win6.png).
  
## Практическая часть
**Порядок компиляции программы:**  
>make
  
В файле lab6.cpp функция slidingWindow реализует алгоритм скользящего окна, функция squareInterp выполняет квадратичную интерполяцию, функция coordinateDescent выполняет покоординатный спуск.  
В файле extraMethods.cpp представлены функции, реализующие дополнительные методы поиска экстремума функции: randomSearch выполняет метод случайного поиска, NelderMead - метод Нелдера-Мида, Powell - метод Пауэлла, HookeJeeves - метод Хука-Дживса, Rosenbrock - метод Розенброка.  
В файле auxiliaryFunctions.cpp находятся вспомогательные функции: func вычисляет значение заданной функции в точке, norm вычисляет норму вектора, normalize нормализирует заданный вектор, sum вычисляет сумму двух векторов, первая функция proizv умножает вектор на число, вторая - вычисляет скалярное произведение двух векторов, функции argmin и argmax выбирают из заданных аргументов тот, при котором функция имеет наименьшее и наибольшее значение соответственно.  
  
**Вывод программы:**  
>Coordinate descent method  
>9 iterations  
>Maximum is on point (27.8, 15.8001)  
>fmax = 294.533  
>  
>Random search method (27.8043, 15.8027)  
>Nelder - Mead method (27.8, 15.8)  
>Powell method (27.8, 15.8)  
>Hooke - Jeeves method (27.8, 15.8)  
>Rosenbrock method (27.8002, 15.7995)  
  
### Результаты
В результате работы программы у функции  
![eq](equations/func2.png)  
был найден экстремум в точке (27.8, 15.8001) (начальная точка (9, 25)) за 9 итераций методом покоординатного спуска с точностью ![eq](equations/eq1.png). Ниже приведен рисунок с изображением линий уровня анализируемой функции и траектория поиска экстремума:  
![img](equations/graph1.png)  
Также были реализованы некоторые дополнительные методы поиска экстремума функции.  
Метод случайного поиска дал результат (27.8043, 15.8027)  
Метод Нелдера-Мида дал результат (27.8, 15.8)  
Метод Пауэлла дал результат (27.8, 15.8)  
Метод Хука-Дживса дал результат (27.8, 15.8)  
Метод Розенброка дал результат (27.8002, 15.7995)  
