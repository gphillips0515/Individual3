# Grace Phillips
# Mr. Messner
# October 18th, 2021


# imports
import os
import sys
from Bio import pairwise2
from Bio.Seq import Seq
from Bio.pairwise2 import format_alignment
import linecache

# collecting user input
userInput = input("Please enter any word: ")
whitespace = " "

# This will throw an error message if the word entered has a space or if the user tries to input more than 1 word and
# restart the program.
if whitespace in userInput:
    print("This can only be one word, please try again. ")
    os.execl(sys.executable, sys.executable, *sys.argv)

# Outputs what the user entered as input
else:
    print('Your word is: ', userInput, '\n')

# Reads in the file of 6 words
temp = open("randomWords.txt", 'r').read().split('\n')

# Storing each word in its own variable
firstWord = linecache.getline('randomwords.txt', 1)
secondWord = linecache.getline('randomwords.txt', 2)
thirdWord = linecache.getline('randomwords.txt', 3)
fourthWord = linecache.getline('randomwords.txt', 4)
fifthWord = linecache.getline('randomwords.txt', 5)
sixthWord = linecache.getline('randomwords.txt', 6)




# This next section will compare the inputted word from the user and compare it with each of the words from the file
# Using the pairwise2 library, we are able to compare the words and see what they have in common, along with scoring
# how similar the words are to each other.


# This is how points will be determined:
# Identical characters are given 2 points, 1 point is deducted for each non-identical character.
# 0.5 points are deducted when opening a gap, and 0.1 points are deducted when extending it.

print("COMPARISON 1:\n")
userAnswer = Seq(userInput)
firstLine = Seq(firstWord)
for firstAlignment in pairwise2.align.globalms(userAnswer, firstLine, 2, -1, -.5, -.1):
    print(format_alignment(*firstAlignment))

print("COMPARISON 2:\n")
userAnswer = Seq(userInput)
secondLine = Seq(secondWord)
for secondAlignment in pairwise2.align.globalms(userAnswer, secondLine, 2, -1, -.5, -.1):
    print(format_alignment(*secondAlignment))

print("COMPARISON 3: \n")
userAnswer = Seq(userInput)
thirdLine = Seq(thirdWord)
for thirdAlignment in pairwise2.align.globalms(userAnswer, thirdLine, 2, -1, -.5, -.1):
    print(format_alignment(*thirdAlignment))

print("COMPARISON 4: \n")
userAnswer = Seq(userInput)
fourthLine = Seq(fourthWord)
for fourthAlignment in pairwise2.align.globalms(userAnswer, fourthLine, 2, -1, -.5, -.1):
    print(format_alignment(*fourthAlignment))

print("COMPARISON 5:\n")
userAnswer = Seq(userInput)
fifthLine = Seq(fifthWord)
for fifthAlignment in pairwise2.align.globalms(userAnswer, fifthLine, 2, -1, -.5, -.1):
    print(format_alignment(*fifthAlignment))

print("COMPARISON 6:\n")
userAnswer = Seq(userInput)
sixthLine = Seq(sixthWord)
for sixthAlignment in pairwise2.align.globalms(userAnswer, sixthLine, 2, -1, -.5, -.1):
    print(format_alignment(*sixthAlignment))

# This takes the scores of each Alignment and stores them within an array
ranking = [firstAlignment.score, secondAlignment.score, thirdAlignment.score,
           fourthAlignment.score, fifthAlignment.score, sixthAlignment.score]

# Sorting the scores
ranking.sort()

# Making the scores show from greatest to least
ranking.reverse()

print("Here are the scores ordered from the most similar to the least similar of your chosen word:", userInput, "\n")
print(ranking)
