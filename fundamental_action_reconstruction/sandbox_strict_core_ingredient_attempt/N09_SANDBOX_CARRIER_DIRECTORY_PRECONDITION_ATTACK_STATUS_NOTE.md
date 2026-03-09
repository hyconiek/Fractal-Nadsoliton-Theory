# N09 Sandbox Carrier-Directory Precondition Attack Status Note

Status: `N09_SANDBOX_CARRIER_DIRECTORY_PRECONDITION_ATTACK_STATUS_NOTE_NO_FALSE_PASS`
As of: `2026-03-09`

## Strongest honest conclusion

The sandbox attack succeeds in one narrow but real sense:

```text
the carrier-directory precondition is no longer a live blocker for the
provider-carrier lane
```

This is stronger than `N08` because the live gap is now one step lower.
The sandbox now has:

1. pair-indexed provider-carrier candidate,
2. filename/path convention,
3. minimal template content,
4. actual carrier directory,
5. but still no provider file and no provider emission.

## What did not happen

This attack did **not** derive:

1. actual provider file creation,
2. actual provider emission,
3. actual `theta_1`, `theta_2`,
4. actual populated `u_1`, `u_2`,
5. actual internal orientation datum,
6. actual `E_orient`,
7. admissible `S_sel_int`,
8. strict-core selector closure.

## Why this attack is still useful

It narrows the live gap one step further:

1. the provider-carrier lane is no longer blocked by missing directory
   infrastructure,
2. the next live blocker is now sharper:
   no provider file has been created yet, and therefore no provider can be
   emitted.

## Honest next move

If the sandbox is continued, the next clean move is:

1. either attack the provider-file creation layer directly,
2. or write the incompatibility boundary that even with directory support the
   provider lane remains non-entering on current sandbox state.
