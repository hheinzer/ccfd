/** \file
 *
 * \brief Provides functions for reading data from the `.ini` parameter file
 *
 * \author hhh
 * \date Sat 21 Mar 2020 10:46:32 AM CET
 */

typedef struct cmd_t cmd_t;

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <stdbool.h>

#include "main.h"

/**
 * \brief A structure used to store the commands, read in from the parameter file
 */
struct cmd_t {
	char key[STRLEN];		/**< the key word of the command */
	char value[STRLEN];		/**< the vale of the command */
	cmd_t *next;			/**< next command */
	cmd_t *prev;			/**< previous command */
};

cmd_t *firstCmd;			/**< first command of the list */

/** \brief Read parameter file and create commands list
 *
 * Read `.ini` file and parse each line into a `cmd_t` object. All `cmd_t`
 * objects are connected in a list of commands starting with `firstCmd`.
 *
 * \param[in] iniFileName The name of the parameter file
 */
void fillCmds(char iniFileName[STRLEN])
{
	FILE *iniFile;
	printf("\nReading Parameter File: '%s'\n", iniFileName);
	iniFile = fopen(iniFileName, "r");
	if (!iniFile) {
		printf("Abort! Could not find specified file.\n");
		exit(1);
	}

	int i, j;
	char line[2 * STRLEN], tmp[STRLEN];
	cmd_t *currCmd = NULL;
	while (fgets(line, sizeof(line), iniFile)) {
		/* read key and make it lowercase */
		j = 0;
		for (i = 0; i < strlen(line); ++i) {
			if (isalnum(line[i])  || (line[i] == '_')) {
				tmp[j++] = tolower(line[i]);
			} else if ((line[i] == '!') || (line[i] == '#') ||
				   (line[i] == '=')) {
				break;
			}
		}
		tmp[j] = '\0';

		/* continue and read value if key was read */
		if (strlen(tmp) > 0) {
			cmd_t *aCmd = malloc(sizeof(cmd_t));
			if (!aCmd) {
				printf("| ERROR: could not allocate command\n");
				exit(1);
			}

			/* store the key */
			strcpy(aCmd->key, tmp);
			aCmd->next = NULL;

			/* read value */
			j = 0;
			for (; i < strlen(line); ++i) {
				if (isalnum(line[i]) || (line[i] == '_') ||
				    (line[i] == '.') || (line[i] == ',') ||
				    (line[i] == '-')) {
					tmp[j++] = line[i];
				} else if ((line[i] == '!') || (line[i] == '#')) {
					break;
				}
			}
			tmp[j] = '\0';

			/* store value */
			strcpy(aCmd->value, tmp);

			/* complete the chain links */
			if (firstCmd == NULL) {
				currCmd = firstCmd = aCmd;
				currCmd->prev = NULL;
			} else {
				currCmd->next = aCmd;
				currCmd->next->prev = currCmd;
				currCmd = currCmd->next;
			}
		}
	}

	fclose(iniFile);
}

/** \brief Delete a single node of the command list.
 *
 * Before deleting the command, the previous command is connected to the next
 * command and vice versa.
 *
 * \param[in] aCmd A pointer to the command that is to be deleted
 */
void deleteCmd(cmd_t *aCmd)
{
	/* link up the previous command to the next, if it exists. Otherwise
	 * we are at the beginning of the list and need to move the first
	 * command one over. */
	if (aCmd->prev) {
		aCmd->prev->next = aCmd->next;
	} else {
		firstCmd = aCmd->next;
	}

	/* link up the next command to the previous command, if it exists */
	if (aCmd->next) {
		aCmd->next->prev = aCmd->prev;
	}

	/* free the memory of the command */
	free(aCmd);
}

/** \brief Find a command in the commands list
 *
 * Find the provided key in the list of commands, and return the address of
 * the corresponding value string. Return NULL if key was not found. Once a
 * key was read from the commands list, it is deleted from the list.
 *
 * \param[in] key Key string of the command to be found
 * \param[out] defMsg String, that indicates if an actual value or the
 *	proposal was returned
 * \param[in] proposal The default value that is used if the key was not
 *	specified
 * \return Pointer to the `value` string, or `NULL`
 */
char *findCmd(const char *key, char defMsg[8], const char *proposal)
{
	/* create a lowercase copy of the key */
	char keyLower[strlen(key) + 1];

	strcpy(keyLower, key);
	for (int i = 0; i < strlen(keyLower); ++i) {
		keyLower[i] = tolower(keyLower[i]);
	}

	/* search for the first matching key in the list */
	char *value = NULL;
	cmd_t *currCmd;
	for (currCmd = firstCmd; currCmd; currCmd = currCmd->next) {
		if (!strcmp(keyLower, currCmd->key)) {
			/* the keys are the same, create enough space for the
			 * value to be stored and break out of the loop */
			value = malloc(strlen(currCmd->value) + 1);
			if (!value) {
				printf("| ERROR: could not allocate value\n");
				exit(1);
			}

			strcpy(value, currCmd->value);
			strcpy(defMsg, "*CUSTOM");
			deleteCmd(currCmd);
			break;
		}
	}

	/* if a proposal is not given, then the key is not optional and is
	 * missing in the .ini file */
	if (!value) {
		if (!proposal) {
			printf("Inifile missing necessary keyword: %s\n",
					keyLower);
			exit(1);
		} else {
			value = malloc(strlen(proposal) + 1);
			if (!value) {
				printf("| ERROR: could not allocate value\n");
				exit(1);
			}

			strcpy(value, proposal);
			strcpy(defMsg, "DEFAULT");
		}
	}

	return value;
}

/** \brief Get a string from the commands list
 *
 * Find the key in the command list and return the corresponding value. If the
 * key is not specified, the proposal will be returned. If the proposal is
 * NULL, but the key is not in the list, an error will be thrown.
 *
 * \param[in] key Key string of the command to be found
 * \param[in] proposal The default value that is used if the key was not
 * \return Pointer to the `value` string, containing the parameter, or the
 *	default value, if the parameter was not specified
 */
char *getStr(const char *key, const char *proposal)
{
	char defMsg[8];
	char *value = findCmd(key, defMsg, proposal);
	printf("| %19s = %27s (%s)\n", key, value, defMsg);
	return value;
}

/** \brief Count how often a key appears
 *
 * Count all occurrences of key in parameter file and return them. If the
 * key is not specified, the proposal will be returned. If the proposal is -1,
 * but the key is not in the list, an error will be thrown.
 *
 * \param[in] key Key string of the command to be found
 * \param[in] proposal The default value that is used if the key was not
 * \return How often the `key` appeared in the command list
 */
int countKeys(const char *key, const int proposal)
{
	/* make lowercase copy of the key */
	char keyLower[strlen(key) + 1];

	strcpy(keyLower, key);
	for (int i = 0; i < strlen(keyLower); ++i) {
		keyLower[i] = tolower(keyLower[i]);
	}

	/* search for occourances of key */
	int count = 0;
	for (cmd_t *currCmd = firstCmd; currCmd; currCmd = currCmd->next) {
		if (!strcmp(keyLower, currCmd->key)) {
			++count;
		}
	}

	if (count == 0) {
		if (proposal == -1) {
			printf("Inifile missing necessary keyword: %s\n",
					keyLower);
			exit(1);
		} else {
			return proposal;
		}
	} else {
		return count;
	}
}

/** \brief Get an integer from the commands list
 *
 * Find the key in the command list and return the corresponding value. If the
 * key is not specified, the proposal will be returned. If the proposal is
 * NULL, but the key is not in the list, an error will be thrown.
 *
 * \param[in] key Key string of the command to be found
 * \param[in] proposal The default value that is used if the key was not
 * \return The value of the parameter, or the default value
 */
int getInt(const char *key, const char *proposal)
{
	char defMsg[8];
	char *valueStr = findCmd(key, defMsg, proposal);
	int value = strtol(valueStr, NULL, 10);
	printf("| %19s = %27d (%s)\n", key, value, defMsg);
	free(valueStr);
	return value;
}

/** \brief Get a double from the commands list
 *
 * Find the key in the command list and return the corresponding value. If the
 * key is not specified, the proposal will be returned. If the proposal is
 * NULL, but the key is not in the list, an error will be thrown.
 *
 * \param[in] key Key string of the command to be found
 * \param[in] proposal The default value that is used if the key was not
 * \return The value of the parameter, or the default value
 */
double getDbl(const char *key, const char *proposal)
{
	char defMsg[8];
	char *valueStr = findCmd(key, defMsg, proposal);
	double value = strtod(valueStr, NULL);
	printf("| %19s = %27g (%s)\n", key, value, defMsg);
	free(valueStr);
	return value;
}

/** \brief Get a boolean from the commands list
 *
 * Find the key in the command list and return the corresponding value. If the
 * key is not specified, the proposal will be returned. If the proposal is
 * NULL, but the key is not in the list, an error will be thrown. The value in
 * the parameter file is accepted as true, if it is a 'T', 't', 'True', or
 * 'true', otherwise it is false.
 *
 * \param[in] key Key string of the command to be found
 * \param[in] proposal The default value that is used if the key was not
 * \return The value of the parameter, or the default value
 */
bool getBool(const char *key, const char *proposal)
{
	char defMsg[8];
	char *valueStr = findCmd(key, defMsg, proposal);
	bool value;
	if (!strcmp(valueStr, "T") || !strcmp(valueStr, "True") ||
		!strcmp(valueStr, "t") || !strcmp(valueStr, "true")) {
		value = true;
	} else {
		value = false;
	}
	printf("| %19s = %27s (%s)\n", key, (value ? "true" : "false"), defMsg);
	free(valueStr);
	return value;
}

/** \brief Get an integer array from the commands list
 *
 * Find the key in the command list and return the corresponding integer array.
 * If the key is not specified, the proposal will be returned. If the proposal
 * is NULL, but the key is not in the list, an error will be thrown.
 *
 * \param[in] key Key string of the command to be found
 * \param[in] N Length of the array that is to be read in
 * \param[in] proposal The default value that is used if the key was not
 * \return A pointer to the value array of the parameter, or the default value
 *	array
 */
int *getIntArray(const char *key, const int N, const char *proposal)
{
	char defMsg[8];
	char *valueStr = findCmd(key, defMsg, proposal);
	int *value = malloc(N * sizeof(int));
	if (!value) {
		printf("| ERROR: could not allocate value\n");
		exit(1);
	}

	char *tok = strtok(valueStr, ",");
	for (int i = 0; i < N; ++i) {
		value[i] = strtol(tok, NULL, 10);
		tok = strtok(NULL, ",");
	}

	char line[STRLEN] = "";
	sprintf(line, "{%d", value[0]);
	for (int i = 1; i < N; ++i) {
		sprintf(line + strlen(line), ", %d", value[i]);
	}
	sprintf(line + strlen(line), "}");

	printf("| %16s[%d] = %27s (%s)\n", key, N, line, defMsg);
	free(valueStr);
	return value;
}

/** \brief Get a double array from the commands list
 *
 * Find the key in the command list and return the corresponding integer array.
 * If the key is not specified, the proposal will be returned. If the proposal
 * is NULL, but the key is not in the list, an error will be thrown.
 *
 * \param[in] key Key string of the command to be found
 * \param[in] N Length of the array that is to be read in
 * \param[in] proposal The default value that is used if the key was not
 * \return A pointer to the value array of the parameter, or the default value
 *	array
 */
double *getDblArray(const char *key, const int N, const char *proposal)
{
	char defMsg[8];
	char *valueStr = findCmd(key, defMsg, proposal);
	double *value = malloc(N * sizeof(double));
	if (!value) {
		printf("| ERROR: could not allocate value\n");
		exit(1);
	}

	char *tok = strtok(valueStr, ",");
	for (int i = 0; i < N; ++i) {
		value[i] = strtod(tok, NULL);
		tok = strtok(NULL, ",");
	}

	char line[STRLEN] = "";
	sprintf(line, "{%g", value[0]);
	for (int i = 1; i < N; ++i) {
		sprintf(line + strlen(line), ", %g", value[i]);
	}
	sprintf(line + strlen(line), "}");

	printf("| %16s[%d] = %27s (%s)\n", key, N, line, defMsg);
	free(valueStr);
	return value;
}

/**
 * \brief Delete all commands in the commands list
 */
void freeCmds(void)
{
	cmd_t *currCmd = firstCmd;
	while (currCmd) {
		if (currCmd->next) {
			currCmd = currCmd->next;
			free(currCmd->prev);
		} else {
			free(currCmd);
			break;
		}
	}
}

/**
 * \brief Print out all remaining commands in the list.
 */
void ignoredCmds(void)
{
	printf("\nThe Following Commands were Ignored:\n");
	for (cmd_t *currCmd = firstCmd; currCmd; currCmd = currCmd->next) {
		printf("| \"%s\" = \"%s\"\n", currCmd->key, currCmd->value);
	}

	/* Remove the remaining commands from memory. */
	freeCmds();
}
