//============================================================================
// Name        : TresEnRaya_Learner.cpp
// Author      : Carlos García Martínez at University of Córdoba (cgarcia@uco.es)
// Version     :
// Copyright   : Licencia Creative Commons Atribución 4.0 Internacional (http://creativecommons.org/licenses/by/4.0/)
//============================================================================

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <utility>
#include <cstdlib>
#include <ctime>
#include <limits>
#include <cmath>
#include <queue>
using namespace std;

//Some initial constants
#define PREF_P1_WIN 2.
#define PREF_P1_TIE 1.
#define PREF_P1_LOSS 2.
#define PREF_P2_WIN 2.
#define PREF_P2_TIE 1.
#define PREF_P2_LOSS 2.
#define MAXOCCURRENCES 100

/**
 * This class stores the statistics about a decision taken on a specific board state
 */
class Decision {

public:
	int**board = NULL; //The board state. This will be a 3x3 matrix with 0 if no tile, 2 for the first player, and 1 for the second player
	int x; //Row (0,1,2) of the decision where to put the next tile
	int y; //Column (0,1,2) of the decision where to put the next tile

	//Statistics
	int numTies;
	int numWins;
	int numLost;

public:

	/**
	 * Copy constructor
	 */
	Decision(Decision &other) {
		x = other.x;
		y = other.y;
		numTies = other.numTies;
		numWins = other.numWins;
		numLost = other.numLost;
		board = new int*[3];

		for (int i = 0; i < 3; i++) {
			board[i] = new int[3];
			for (int j = 0; j < 3; j++) {
				board[i][j] = other.board[i][j];
			}
		}
	}

	/**
	 * Constructor
	 */
	Decision() {
		x = 0;
		y = 0;
		numTies = 0;
		numWins = 0;
		numLost = 0;
		board = new int*[3];

		for (int i = 0; i < 3; i++) {
			board[i] = new int[3];
			for (int j = 0; j < 3; j++) {
				board[i][j] = 0;
			}
		}
	}

	/**
	 * Constructor that reads the information from a stream
	 */
	Decision(istream &stream) {
		board = new int*[3];

		for (int i = 0; i < 3; i++)
			board[i] = new int[3];

		read(stream);
	}

	/**
	 * Basic constructor with arguments
	 */
	Decision(int **board, int x, int y) {
		this->x = x;
		this->y = y;
		this->board = new int*[3];

		for (int i = 0; i < 3; i++) {
			this->board[i] = new int[3];
			for (int j = 0; j < 3; j++) {
				this->board[i][j] = board[i][j];
			}
		}

		numTies = 0;
		numWins = 0;
		numLost = 0;
	}

	/**
	 * Destructor
	 */
	~Decision() {

		for (int i = 0; i < 3; i++) {
			delete[] board[i];
		}

		delete board;
	}

	/**
	 * This function writes the object onto a stream
	 */
	void write(ostream &stream) {

		for (int i = 0; i < 3; i++) {
			for (int j = 0; j < 3; j++) {
				stream << board[i][j] << " ";
			}
		}

		stream << x << " " << y << " " << numTies << " " << numWins << " "
				<< numLost << endl;
	}

	/**
	 * This function initialise the object's variables from a stream
	 */
	void read(istream &stream) {

		for (int i = 0; i < 3; i++) {
			for (int j = 0; j < 3; j++) {
				stream >> board[i][j];
			}
		}

		stream >> x;
		stream >> y;
		stream >> numTies;
		stream >> numWins;
		stream >> numLost;
	}

	/**
	 * This function compares if the object represents the same state as that passed by argument (board and decision)
	 */
	bool matches(int **board, int x, int y) {

		if (this->x != x || this->y != y)
			return false;

		for (int i = 0; i < 3; i++) {
			for (int j = 0; j < 3; j++) {
				if (this->board[i][j] != board[i][j])
					return false;
			}
		}

		return true;
	}

	/**
	 * This function compares if the object considers the same board state as that passed by argument
	 */
	bool matchesBoard(int **board) {

		for (int i = 0; i < 3; i++) {
			for (int j = 0; j < 3; j++) {
				if (this->board[i][j] != board[i][j])
					return false;
			}
		}

		return true;
	}

};

/**
 * This function reads a set of Decision objects from a stream
 */
void readMemory(istream &stream, vector<pair<Decision*, bool>> &memory) {

	while (stream.peek() == ' ')
		stream.get();

	while (stream.peek() != EOF) {
		memory.push_back(make_pair(new Decision(stream), false));

		while (stream.peek() == ' ' && stream.peek() != EOF)
			stream.get();
	}

	if (memory.size() > 0){
		Decision *dec = memory.back().first;
		delete dec;
		memory.pop_back();
	}
}

/**
 * This function writes a set of Decision objects onto a stream
 */
void writeMemory(ostream &stream, vector<pair<Decision*, bool>> &memory) {

	for (auto pair : memory) {
		pair.first->write(stream);
	}
}

/**
 * This function prints a board state to the standard output
 */
void printBoard(int **board) {
	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++) {
			cout << board[i][j] << " ";
		}

		cout << endl;
	}
}

/**
 * This function returns a set of posible decisions to be taken from a board state
 */
void getAlternatives(int **board, vector<pair<int, int>> &alternatives) {

	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++) {
			if (board[i][j] == 0) {
				alternatives.push_back(make_pair(i, j));
			}
		}
	}
}

/**
 * This function looks for a decision (board + decision) in a set of decisions
 */
int findState(int **board, int x, int y,
		vector<pair<Decision*, bool>> &memory) {

	int index = 0;

	for (auto aState : memory) {
		if (aState.first->matches(board, x, y)) {
			return index;
		}
		index++;
	}

	return -1;
}

/**
 * This function randomly selects a number in [0, probabilities.size() - 1] according to the values in that vector. Larger values are selected with a higher probability
 */
int chance(vector<double>&probabilities) {

	double sum = 0;

	for (auto prob : probabilities)
		sum += prob;

	double sample = (((double) rand()) / RAND_MAX) * sum;
	int selected = 0;
	int probSize = probabilities.size();

	while (sample >= 0 && selected < probSize) {
		sample -= probabilities.at(selected);
		selected++;
	}

	selected--;
	return selected;
}

/**
 * This function determines if the game has ended
 */
int testWinState(int **board) {

	int winner = 0;

	//Row checks
	for (int i = 0; i < 3 && winner == 0; i++) {
		if (board[i][0] != 0) {
			winner = board[i][0];
			for (int j = 1; j < 3 && winner != 0; j++) {

				if (board[i][j] != winner)
					winner = 0;
			}
		}
	}

	//Column checks
	if (winner == 0) {
		for (int j = 0; j < 3 && winner == 0; j++) {
			if (board[0][j] != 0) {
				winner = board[0][j];
				for (int i = 0; i < 3 && winner != 0; i++) {

					if (board[i][j] != winner)
						winner = 0;
				}
			}
		}
	}

	//First diagonal check
	if (winner == 0) {
		if (board[0][0] != 0) {
			winner = board[0][0];
			for (int k = 1; k < 3 && winner != 0; k++) {

				if (board[k][k] != winner)
					winner = 0;
			}
		}
	}

	//Second diagonal check
	if (winner == 0) {
		if (board[2][0] != 0) {
			winner = board[2][0];
			for (int k = 1; k < 3 && winner != 0; k++) {

				if (board[2 - k][k] != winner)
					winner = 0;
			}
		}
	}

	return winner;
}

/**
 * Main function
 */
int main(int argc, char *argv[0]) {

	//Random numbers initialisation
	srand(time(0));
	time_t initTime;
	time_t currentTime;

#ifdef __DEBUG__
	Decision dec1;
	dec1.board[0][0] = 0;
	dec1.board[0][1] = 0;
	dec1.board[0][2] = 2;
	dec1.board[1][0] = 1;
	dec1.board[1][1] = 2;
	dec1.board[1][2] = 1;
	dec1.board[2][0] = 2;
	dec1.board[2][1] = 1;
	dec1.board[2][2] = 2;
	testWinState(dec1.board);
	dec1.x = 0;
	dec1.y = 2;
#endif

	//
	queue<int> result;

	if (argc != 5) {
		cout << argv[0]
				<< " <memfile> <seconds> <probSmoother_P1> <probSmoother_P2>"
				<< endl;
		cout
				<< "<probSmoother_PX> is the probability smoother of the first player. It should be between 0 and 100, with 0 for no smoothing (strict) and 100 for too much smoothing (random)"
				<< endl;
		return 0;
	}

	int max_seconds = 1;
	int p1_Wins = 0;
	int numTies = 0;
	int impresedComputWins = -1;
	int impressedPseudoHumanWins = -1;
	int impressedTies = -1;
	int p2_Wins = 0;

	max_seconds = atoi(argv[2]);
	double p1_smoothing = atof(argv[3]);
	double p2_smoothing = atof(argv[4]);

	if (p1_smoothing > 100 || p1_smoothing < 0) {
		cerr << "P1 smoothing set to 0";
		p1_smoothing = 0;
	}

	if (p2_smoothing > 100 || p2_smoothing < 0) {
		cerr << "P2 smoothing set to 0";
		p2_smoothing = 0;
	}

	double ksmoothing_p1 = p1_smoothing;
	double k3smoothing_p1 = 3. * p1_smoothing;
	double ksmoothing_p2 = p2_smoothing;
	double k3smoothing_p2 = 3. * p2_smoothing;

	//Read the memory from the file, if exists
	vector<pair<Decision*, bool>> memory;
	ifstream file(argv[1], ifstream::in);
	if (!file.is_open()) {
		cerr << "Creating file " << argv[1] << endl;
	} else {
		readMemory(file, memory);
		file.close();
	}

	//The board
	int **board = new int*[3];
	for (int i = 0; i < 3; i++)
		board[i] = new int[3];

	time(&initTime);
	time(&currentTime);
	long maxSeconds = max_seconds;
	unsigned long iIterations = 0;

	//Main loop
	while (difftime(currentTime, initTime) < maxSeconds) {

		int sizeMemory = memory.size();

		//reset the chosen decision from the memory
		for (int i = 0; i < sizeMemory; i++)
			memory.at(i).second = false;

		//boolean variables to detect if the first player won, lost, or tied
		bool iLost = false;
		bool iWon = false;
		bool thereAreAlternatives = true;

		//reset the board
		for (int i = 0; i < 3; i++) {
			for (int j = 0; j < 3; j++) {
				board[i][j] = 0;
			}
		}

		//While the game had not ended
		while (iLost == false && iWon == false && thereAreAlternatives) {
			int x, y;

			//**********************
			//First player turn

			//Get the alternatives
			vector<pair<int, int>> alternatives;
			vector<double> probabilities;
			getAlternatives(board, alternatives);

			//Tie condition
			if (alternatives.size() <= 0) {
				thereAreAlternatives = false;
				continue;
			}

			//Compute the quality of each alternative
			for (auto alternative : alternatives) {

				//Retrive the statistics for each alternative
				int index = findState(board, alternative.first,
						alternative.second, memory);
				int numWins = 0;
				int numTies = 0;
				int occurrences = 0;
				int numLost = 0;

				if (index >= 0) {
					Decision *state = memory.at(index).first;
					numWins = state->numWins;
					numTies = state->numTies;
					numLost = state->numLost;
					occurrences = numWins + numLost + numTies;
				} else {
					//If the alternative is not in the memory, create and initialise it, and push it into the memory
					Decision *dec = new Decision(board, alternative.first,
							alternative.second);
					dec->numTies = MAXOCCURRENCES / 3;
					dec->numWins = MAXOCCURRENCES / 3;
					dec->numLost = MAXOCCURRENCES / 3;
					memory.push_back(make_pair(dec, false));
				}

				//Compute the probabilities
				double probWin = (numWins + ksmoothing_p1)
						/ (occurrences + k3smoothing_p1);
				double probLoss = (numLost + ksmoothing_p1)
						/ (occurrences + k3smoothing_p1);
				double probTie = (numTies + ksmoothing_p1)
						/ (occurrences + k3smoothing_p1);

				//Compute the quality of the alternative as the probability of winning, powered by a factor, times one less the probability of lossing, also powered by a factor
				double quality = pow(probWin, PREF_P1_WIN)
						* pow((1 - probLoss), PREF_P1_LOSS);
				probabilities.push_back(quality);
			}

			//Choose one alternative
			int election = chance(probabilities);
			election = (election >= probabilities.size()) ? 0 : election;
			x = alternatives.at(election).first;
			y = alternatives.at(election).second;
			int index = findState(board, x, y, memory);
			memory.at(index).second = true;

			//Put the tile
			board[x][y] = 2;

#ifdef __DEBUG__
			printBoard(board);
#endif

#ifdef __DEBUG__
			if (dec1.matchesBoard(board)) {
				cout << "hola" << endl;
			}
#endif

			//Test if game ended
			if (testWinState(board) == 2) {
				p1_Wins++;
				iWon = true;
				continue;
			}

			//Check the tie condition
			alternatives.clear();
			probabilities.clear();
			getAlternatives(board, alternatives);

			if (alternatives.size() <= 0) {
				thereAreAlternatives = false;
				continue;
			}

			//*****************************
			//P2'S TURN

			//Compute the quality of each alternative
			for (auto alternative : alternatives) {

				//Retrieve the statistics of each alternative
				int index = findState(board, alternative.first,
						alternative.second, memory);
				int numWins = 0;
				int numTies = 0;
				int occurrences = MAXOCCURRENCES;
				int numLost = 0;

				if (index >= 0) {
					Decision *state = memory.at(index).first;
					numWins = state->numLost; //This player uses the statistics the other way round!!
					numTies = state->numTies;
					numLost = state->numWins;
					occurrences = numWins + numLost + numTies;
				} else {
					//If the alternative is not in the memory, create and initialise it, and push it into the memory
					Decision *dec = new Decision(board, alternative.first,
							alternative.second);
					dec->numTies = MAXOCCURRENCES / 3;
					dec->numWins = MAXOCCURRENCES / 3;
					dec->numLost = MAXOCCURRENCES / 3;
					memory.push_back(make_pair(dec, false));
				}

				//Compute the probabilities
				double probWin = (numWins + ksmoothing_p2)
						/ (occurrences + k3smoothing_p2);
				double probLoss = (numLost + ksmoothing_p2)
						/ (occurrences + k3smoothing_p2);
				double probTie = (numTies + ksmoothing_p2)
						/ (occurrences + k3smoothing_p2);

				//Compute the quality of the alternative
				double quality = pow(probWin, PREF_P1_WIN)
						* pow((1 - probLoss), PREF_P1_LOSS);
				probabilities.push_back(quality);
			}

			//Choose one alternative and mark it as used
			election = chance(probabilities);
			election = (election >= probabilities.size()) ? 0 : election;
			x = alternatives.at(election).first;
			y = alternatives.at(election).second;
			index = findState(board, x, y, memory);
//					x++;
//					y++;
			memory.at(index).second = true;

//			x--;
//			y--;

			//Put the tile
			board[x][y] = 1;

#ifdef __DEBUG__
			printBoard(board);
#endif

#ifdef __DEBUG__
			if (dec1.matchesBoard(board)) {
				cout << "hola" << endl;
			}
#endif

			//Test the end of the game
			if (testWinState(board) == 1) {
				p2_Wins++;
				iLost = true;
				continue;
			}
		}

		//Update the statistics in the memory once the game ended
		for (auto pair : memory) {
			if (pair.second) {
				int total = pair.first->numWins;
				total += pair.first->numTies;
				total += pair.first->numLost;
				total++;

				if (iWon) {
					if (pair.first->numWins < MAXOCCURRENCES)
						pair.first->numWins++;
					if (total > MAXOCCURRENCES) {
						if (pair.first->numTies > 0)
							pair.first->numTies--;
						else if (pair.first->numLost > 0)
							pair.first->numLost--;
					}
				} else if (iLost) {
					if (pair.first->numLost < MAXOCCURRENCES)
						pair.first->numLost++;
					if (total > MAXOCCURRENCES) {
						if (pair.first->numTies > 0)
							pair.first->numTies--;
						else if (pair.first->numWins > 0)
							pair.first->numWins--;
					}
				} else {
					if (pair.first->numTies < MAXOCCURRENCES)
						pair.first->numTies++;
					if (total > MAXOCCURRENCES) {
						if (pair.first->numWins > 0)
							pair.first->numWins--;
						else if (pair.first->numLost > 0)
							pair.first->numLost--;
					}
				}
			}
		}

		time(&currentTime);

		iIterations++;

		if (iWon) {
			result.push(2);
		} else if (iLost) {
			result.push(1);
		} else {
			numTies++;
			result.push(0);
		}

		if (result.size() > 1000) {
			int oldestResult = result.front();
			result.pop();

			if (oldestResult == 2)
				p1_Wins--;
			else if (oldestResult == 1)
				p2_Wins--;
			else
				numTies--;
		}


		//Print the current statistics
		if ((abs(p1_Wins - impresedComputWins)
				+ abs(p2_Wins - impressedPseudoHumanWins)
				+ abs(numTies - impressedTies)) > 50) {

			cout << iIterations << " " << p1_Wins << " " << p2_Wins
					<< " " << numTies << endl;
			impresedComputWins = p1_Wins;
			impressedPseudoHumanWins = p2_Wins;
			impressedTies = numTies;
//			cout << "Comput: " << computWins << " Human: " << p2_Wins
//					<< " Ties: " << numTies << endl;
//		cout << "Comput: " << ((100 * computWins) / iIterations) << " Human: "
//				<< ((100 * pseudoHumanWins) / iIterations) << " Ties: "
//				<< ((100 * (iIterations - computWins - pseudoHumanWins))
//						/ iIterations) << endl;
		}
	}

	ofstream fileOutput(argv[1], ofstream::out);
	writeMemory(fileOutput, memory);
	fileOutput.close();
	return 0;
}
